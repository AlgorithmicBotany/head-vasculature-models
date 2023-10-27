#include "include_all.h"

#include "lpfg_continuous_panel_utility.h"
#include <lpfgall.h>
// ordering: must be after <lpfgall.h>
#include "lpfg_draw_utility.h"

struct Model {
  virtual ~Model() = default;
  virtual void init() = 0;
  virtual void step() = 0;
  virtual void draw() = 0;
};
std::unique_ptr<Model> model;

namespace phyl = modelling::phyllotaxis;
namespace vasc = vasculature;

struct BellisVasculature : Model {
  using V3 = math::Vec3f;

  using MRRodkaew = vasc::MinimalResistanceRodkaew;
  using Vasculature = MRRodkaew::Vasculature;

  Vasculature veins;

  // vein segmentations and groups
  // warning: no vein nodes are ever deleted, so it is fine to store node_ptr
  // here, however, some 'book keeping' would be necessary if veins (nodes) are
  // ever deleted.
  Vasculature::node_refs primordia;
  Vasculature::node_refs stemBase;
  Vasculature::node_refs bracts;
  Vasculature::node_refs florets;
  Vasculature::node_refs incipientVeins;

  MRRodkaew::DiffuseAuxinPoints diffuseAuxin;

  // Reticulation (refactor maybe?)
  Vasculature::node_refs reticulateNodesToBeProcessed;
  Vasculature::node_refs reticulateNodes;

  Vasculature::node_refs subReticulateNodes;

  struct Params {
    float r_groundTissue{1.75f};
    float r_veins{1.f};
    float veinFullGrownAge{250.f};

    float mergeDistance{0.04f};
    float stepDistance{0.02f};
    float killDistance{0.2f};

    float veinJitterFrequency{0.5f};
    float veinJitterStrength{0.5f};

    float sectorialPropensity{0.1f};

    float floretInitialKick{1.f};
    float bractInitialKick{1.f};

    bool useReticulation{false};
    float reticulationPropensity{0.1f};
    int reticulationAgeTrigger{1};

    float shieldingSolidAngle{0.1f};
  } params;

  std::vector<phyl::TimestampPrimordia> orderedPrims;
  phyl::GrowingReceptacle growingReceptacle;

  int frame;
  int iterPerStep{10};

  struct Viz {
    struct Color {
      int adaxial{6};
      int bract{7};
      int bractPrim{21};
      int stem{7};
      int primordia{8};
      int floret{8};
      int floretPrim{20};
      int auxin{9};
      int backside{56};
    } color;
    struct Size {
      float floret{0.5f};
      float bract{0.5f};
      float stem{0.5f};
      float primordia{0.5f};
    } size;
    struct VeinWidth {
      float floret{1.f};
      float adaxial{1.f};
      float abaxial{1.f};
    } veinWidth;

    struct Render {
      bool oldestNeighbours{true};
      bool closestNeighbours{false};
      bool dashedContour{false};
      bool diffuseAuxin{false};
      enum class Style { Normal, EffectiveResistance, AgeResistance } style;
    } render;

    int neighbourCombo{0};
    float neighbourLength{0.5f};
    float veinScale{1.f};
    float veinExponent{0.35f};
    float tipWidth{1.f};

    int veinAgeFadeIn{100};
    bool continuousRereadPanel{false};

    float backsideDepth{0.01f};

  } viz;

  int colorID(MRRodkaew::NodeType type) {
    using Type = MRRodkaew::NodeType;
    switch (type) {
    case Type::FLORET:
      return viz.color.floret;
    case Type::ABAXIAL:
      return viz.color.bract;
    case Type::ADAXIAL:
      return viz.color.adaxial;
    default:
      return 0;
    }
  }
  int colorID(Vasculature::node_ptr node) {
    if (node->data.position.z < -viz.backsideDepth)
      return viz.color.backside;
    return colorID(node->data.type);
  }

  float veinTypeWidth(MRRodkaew::NodeType type) {
    using Type = MRRodkaew::NodeType;
    switch (type) {
    case Type::FLORET:
      return viz.veinWidth.floret;
    case Type::ABAXIAL:
      return viz.veinWidth.abaxial;
    case Type::ADAXIAL:
      return viz.veinWidth.adaxial;
    default:
      return 1.f;
    }
  }
  float veinTypeWidth(Vasculature::node_ptr node) {
    return veinTypeWidth(node->data.type);
  }
  float veinAgeFadeIn(Vasculature::node_ptr node) {
    // "pop" into existence
    // return node->data.width;
    // "fade" into existance
    return node->data.width *
           func(vein_width_fadein_f, node->data.age / float(viz.veinAgeFadeIn));
  }
  float veinWidth(Vasculature::node_ptr node) { return veinAgeFadeIn(node); }

  //
  // PHYLLOTAXIS
  //
  std::vector<phyl::TimestampPrimordia> gatherNextPrimordia() {
    auto frameA = frame;
    auto frameB = frame + 1;
    return phyl::gatherTimestampPrimordiaBetween(frameA, frameB, orderedPrims);
  }

  bool inAngleFilter(float angle) {
    const float showAngle = 0.5;

    return angle > 0 && angle < M_PI;

    return std::abs(angle - M_PI) < showAngle ||
           std::abs(angle - 2 * M_PI) < showAngle ||
           std::abs(angle) < showAngle;
  }
  void debug_filterDiffuseAuxin(std::vector<V3> &points) {
    for (int i = 0; i < points.size(); ++i) {
      auto const &p = points[i];
      auto angle = std::atan2(-p.z, p.x) + M_PI;
      if (!inAngleFilter(angle))
        utility::swap_pop(points, i--);
    }
  }

  void debug_filterPrimordia(std::vector<phyl::TimestampPrimordia> &primordia) {
    auto iter = std::stable_partition(
        std::begin(primordia), std::end(primordia),
        // C++17 [=](auto const &prim) -> bool {
        [=](phyl::TimestampPrimordia const &prim) -> bool {
          return inAngleFilter(prim.angle);
        });
    primordia.erase(iter, std::end(primordia));
  }

  void debug_build3ClosestAndOldestPrims(Vasculature::node_ptr node,
                                         Vasculature::node_refs prims) {
    std::sort(std::begin(prims), std::end(prims),
              // c++ 17 [=](auto a, auto b) {
              [=](Vasculature::node_ptr a, Vasculature::node_ptr b) {
                auto dA = distance(node->data.position, a->data.position);
                auto dB = distance(node->data.position, b->data.position);
                return dA < dB;
              });

    int id = 0;
    for (; id < 3 && id < prims.size(); ++id) {
      if (distance(node->data.position, prims[id]->data.position) > 0.6f)
        break;
      node->data.closestNeighs.push_back(&(prims[id]->data));
    }

    prims.resize(std::min(id, 3));

    std::sort(std::begin(prims), std::end(prims),
              //   C++ 17   [=](auto a, auto b)
              [=](Vasculature::node_ptr a, Vasculature::node_ptr b) {
                return a->data.age < b->data.age;
              });

    for (int id = 0; id < 3 && id < prims.size(); ++id) {
      node->data.oldestNeighs.push_back(&(prims[id]->data));
    }
  }

  void debug_buildSolidAngleClustering(Vasculature::node_ptr node,
                                       Vasculature::node_refs primorida) {

    auto p = node->data.position;
    std::vector<V3> points;
    points.reserve(primorida.size());

    for (auto prim : primorida) {
      points.push_back(prim->data.position);
    }

    auto prims =
        // clustering::selectSolidAngleClusterAround(p, points, 1.8f, 2.f *
        // M_PI);
        clustering::selectSolidAngleClusterAround(p, points, 0.49f, 2.f * M_PI);

    Vasculature::node_refs considered;
    int id = 0;
    for (; id < prims.size(); ++id) {
      auto pos = prims[id].p * prims[id].w + p;
      // if (distance(p, pos) > 0.6f)
      //  break;
      float closestDist = 99999.f;
      Vasculature::node_ptr closest;
      for (auto other : veins) {
        auto dist = distance(other->data.position, pos);
        if (dist < closestDist) {
          closestDist = dist;
          closest = other;
        }
      }
      considered.push_back(closest);
    }

    auto minEffRes = 0.f;
    Vasculature::node_ptr minNode = nullptr;
    for (auto c : considered) {
      // node->data.oldestNeighs.push_back(&(c->data));

      auto dis = distance(p, c->data.position);
      auto effRes = dis * (params.r_groundTissue + c->data.resistanceAtNode) +
                    c->data.effectiveResistance;
      effRes = c->data.age;
      if (effRes > minEffRes) {
        minEffRes = effRes;
        minNode = c;
      }
    }

    // if (!considered.empty() && minNode != nullptr)
    //  node->data.oldestNeighs.push_back(&(minNode->data));
  }

  void addPrimordia() {

    auto receptaclePrims = gatherNextPrimordia();

    // <dbg>
    // debug_filterPrimordia(receptaclePrims);
    // </dbg>

    for (auto const &prim : receptaclePrims) {
      auto primPoint = phyl::receptaclePoint(prim, growingReceptacle);

      auto node = veins.addNode();
      node->data.position = primPoint.p;
      node->data.direction = primPoint.n;
      node->data.resistanceAtNode = params.r_groundTissue;
      node->data.effectiveResistance = math::maxValue<float>();

      // <dbg
      debug_build3ClosestAndOldestPrims(node, primordia);
      debug_buildSolidAngleClustering(node, veins.nodes());
      // </dgb

      primordia.push_back(node);
      incipientVeins.push_back(node);

      switch (prim.type) {
      case phyl::TimestampPrimordia::Type::FLORET: {
        florets.push_back(node);
        // node->data.direction.y -= 1.f;
        node->data.direction *=
            func(floret_depth_f, growingReceptacle.currentTime()) *
            params.floretInitialKick;
        node->data.type = MRRodkaew::NodeType::FLORET;
        node->data.baseResistance = 5.f;
        node->data.age = 0;
        reticulateNodesToBeProcessed.push_back(node);
        break;
      }
      case phyl::TimestampPrimordia::Type::BRACT: {
        bracts.push_back(node);
        node->data.direction *= params.bractInitialKick;
        node->data.type = MRRodkaew::NodeType::ABAXIAL;
        node->data.baseResistance = 1.f;

        node->data.age = 0;
        reticulateNodesToBeProcessed.push_back(node);

        // make bract laterals
        if (false) {
          auto lateralPrim = prim;
          lateralPrim.angle += 0.1f;
          auto lateralPoint =
              phyl::receptaclePoint(lateralPrim, growingReceptacle);

          auto lateral = veins.addNode();
          lateral->data.position = lateralPoint.p;
          lateral->data.direction = lateralPoint.n;
          lateral->data.resistanceAtNode = params.r_groundTissue;
          lateral->data.effectiveResistance = math::maxValue<float>();

          bracts.push_back(lateral);
          lateral->data.direction *= params.bractInitialKick;
          lateral->data.type = MRRodkaew::NodeType::ABAXIAL;

          lateral->data.age = -55;
          incipientVeins.push_back(lateral);
        }
        // make bract laterals
        if (false) {
          auto lateralPrim = prim;
          lateralPrim.angle -= 0.1f;
          auto lateralPoint =
              phyl::receptaclePoint(lateralPrim, growingReceptacle);

          auto lateral = veins.addNode();
          lateral->data.position = lateralPoint.p;
          lateral->data.direction = lateralPoint.n;
          lateral->data.resistanceAtNode = params.r_groundTissue;
          lateral->data.effectiveResistance = math::maxValue<float>();

          bracts.push_back(lateral);
          lateral->data.direction *= params.bractInitialKick;
          lateral->data.type = MRRodkaew::NodeType::ABAXIAL;

          lateral->data.age = -55;
          incipientVeins.push_back(lateral);
        }

        break;
      }
      }
    }
  }

  void addDiffuseAuxin() {

    int const sampleCount = 20;
    float const packingRadius = 0.1f;

    auto t = growingReceptacle.currentTime();

    auto arcLength_diffuseAuxin = func(diffuse_auxin_f, t);
    auto newDiffuseAuxin = phyl::generateNewAdaxialAuxin(
        sampleCount, packingRadius, diffuseAuxin, arcLength_diffuseAuxin,
        growingReceptacle.currentCage());

    // <dbg>
    debug_filterDiffuseAuxin(newDiffuseAuxin);
    // </dbg>
    utility::append(newDiffuseAuxin, diffuseAuxin);
  }

  Vasculature::node_refs gatherVeinTips() {
    Vasculature::node_refs veinTips;

    using namespace std;

    auto iter = partition(
        begin(incipientVeins), end(incipientVeins),
        [](Vasculature::node_ptr node) { return node->data.age < 0; });

    veinTips.reserve(distance(iter, end(incipientVeins)));
    copy(iter, end(incipientVeins), back_inserter(veinTips));
    incipientVeins.erase(iter, end(incipientVeins));

    return veinTips;
  }

  //
  // RECEPTACLE GROWTH AND DEFORMATION
  //
  void receptacleGrowthAndVeinDeformation() {
    if (growingReceptacle.isFullGrown())
      return;

    geometry::insitu_subdivideAtLength(
        veins,
        // c++17 [](auto const & a, auto const & b) -> MRRodkaew::Node {
        [](MRRodkaew::Node const &a,
           MRRodkaew::Node const &b) -> vasc::MinimalResistanceRodkaew::Node {
          vasc::MinimalResistanceRodkaew::Node mid;
          mid.age = a.age;
          mid.width = 0.5f * (a.width + b.width),
          mid.position = math::mid(a.position, b.position),
          mid.resistanceAtNode = a.resistanceAtNode;
          mid.baseResistance = a.baseResistance;

          mid.type = b.type;
          return mid;
        },
        0.1f);

    // 1. pack a contiguous array of points (saving intervals)
    // 2. project into 2D contour region
    // 3. deform via contour cage
    // 4. put deformed points back into there respective containers

    // 1. pack tmp points for deformation
    std::vector<math::Vec3f> tmpPoints;
    auto intervalVeins = utility::intervalAppend(
        veins.nodes(), tmpPoints,
        // c++ 17[](auto node)
        [](Vasculature::node_ptr node) { return node->data.position; });

    auto intervalAuxinPoints = utility::intervalAppend(diffuseAuxin, tmpPoints);

    // 2., 3. deform tmp points from current to next cage
    auto const &currentCage = growingReceptacle.currentCage();
    growingReceptacle.step();
    auto const &nextCage = growingReceptacle.currentCage();
    modelling::insitu_deformVolumeOfRevolutionPoints(tmpPoints, currentCage,
                                                     nextCage);

    // 4. unpack deformed points
    utility::intervalUpdate(intervalVeins, tmpPoints, veins.nodes(),
                            // c++ 17 [](auto a, auto b)
                            [](Vasculature::node_ptr a, math::Vec3f const &b) {
                              a->data.position = b;
                            });
    utility::intervalCopy(intervalAuxinPoints, tmpPoints, diffuseAuxin);
  }

  //
  // VASCULAR DEVELOPMENT
  //
  void makeStemBase() {
    auto const &currentCage = growingReceptacle.currentCage();
    auto point = currentCage.controlPoints.points()[55];
    auto target = math::vec3f::from(point) - V3{0.025f, 0.f, 0.f};
    // auto target = math::vec3f::from(point) - V3{0.025f, 0.f, 0.f};

    for (int i = 0; i < 52; ++i) {
      float t = (i / 52.f) * 360.f;
      auto p = math::rotateAroundY(target, t);

      auto sink = veins.addNode();

      sink->data.position = p;
      sink->data.age = params.veinFullGrownAge;
      sink->data.id = i + 10000;
      sink->data.type = MRRodkaew::NodeType::ABAXIAL;
      sink->data.baseResistance = 1.f;
      sink->data.distanceToBase = 0.f;
      sink->data.effectiveResistance = 0.f;
      sink->data.resistanceAtNode = 0.f;

      stemBase.push_back(sink);
    }
  }

  //
  // VISUALIZATION
  //
  void recalculateVascularWidths() {
    // reset widths
    MRRodkaew::Vasculature::node_refs tips;
    for (auto node : veins) {
      node->data.width = 0.f;

      if (node->data.type == MRRodkaew::NodeType::FLORET) {
        for (auto child : node->children()) {
          if (child->data.type == MRRodkaew::NodeType::ADAXIAL) {
            node->data.type = MRRodkaew::NodeType::ADAXIAL;
            break;
          }
        }
      }

      if (node->isLeaf())
        tips.push_back(node);
    }

    // accumulate from tips to sinks
    for (auto node : tips) {
      node->data.width += viz.tipWidth;

      while (node->hasParents()) {
        node = node->parents()[0];
        node->data.width += viz.tipWidth;
      }
    }

    for (auto node : subReticulateNodes) {
      if (node->data.closestNeighs.empty())
        continue;

      auto width = node->data.closestNeighs.front()->width;
      while (node->hasParents()) {
        node = node->parents()[0];
        node->data.width += width * 0.5f;
      }
    }

    // recalculate widths
    for (auto node : veins) {
      auto &w = node->data.width;
      w = viz.veinScale * veinTypeWidth(node) * std::pow(w, viz.veinExponent);
    }
  }

  void rereadPanel() {
    viz.continuousRereadPanel = cval("CONTINUOUS_REREAD");

    if (!viz.continuousRereadPanel)
      return;

    viz.veinScale = cval("VEIN_SCALE");
    viz.veinExponent = cval("VEIN_EXPONENT");
    viz.tipWidth = cval("VEIN_TIP_WIDTH");

    viz.size.primordia = cval("PRIM_SIZE");
    viz.size.floret = cval("FLORET_PRIM_SIZE");
    viz.size.bract = cval("BRACT_PRIM_SIZE");
    viz.size.stem = cval("STEM_NODE_SIZE");

    viz.veinWidth.floret = cval("FLORET_VEIN_SCALE");
    viz.veinWidth.adaxial = cval("ADAXIAL_VEIN_SCALE");
    viz.veinWidth.abaxial = cval("ABAXIAL_VEIN_SCALE");

    viz.veinScale = cval("VEIN_SCALE");
    viz.veinExponent = cval("VEIN_EXPONENT");
    viz.veinAgeFadeIn = cval("VEIN_AGE_FADE_IN");

    auto renderNeigh = cval("RENDER_NEIGHBOURS");
    viz.render.closestNeighbours = renderNeigh == 1 ? true : false;
    viz.render.oldestNeighbours = renderNeigh == 2 ? true : false;
    viz.neighbourLength = cval("NEIGHBOUR_LENGTH");
    viz.neighbourCombo = cval("NEIGHBOUR_COMBO");

    auto style = int(cval("RENDER_STYLE"));
    switch (style) {
    case 0:
      viz.render.style = Viz::Render::Style::Normal;
      break;
    case 1:
      viz.render.style = Viz::Render::Style::AgeResistance;
      break;
    case 2:
      viz.render.style = Viz::Render::Style::EffectiveResistance;
      break;
    };
    viz.render.dashedContour = cval("RENDER_DASHED") > 0.f ? true : false;
    viz.render.diffuseAuxin = cval("RENDER_DIFFUSE_AUXIN") > 0.f ? true : false;
    viz.backsideDepth = cval("BACK_SIDE_AT_DEPTH");

    float t = growingReceptacle.currentTime();
    params.r_groundTissue = func(r_ground_f, t);
    params.r_veins = func(r_vein_f, t);
    params.veinFullGrownAge = cval("VEIN_AGE_FULLGROWN");

    params.veinJitterFrequency = cval("VEIN_JITTER_FREQ");
    params.veinJitterStrength = cval("VEIN_JITTER_STRENGTH");

    params.sectorialPropensity = cval("SECTORIAL_PROPENSITY");

    params.floretInitialKick = cval("FLORET_KICK");
    params.bractInitialKick = cval("BRACT_KICK");
    params.useReticulation = bool(cval("USE_RETICULATION"));
    params.reticulationPropensity = cval("RETICULATION_PROPENSITY");
    params.reticulationAgeTrigger = cval("RETICULATION_AGE_TRIGGER");

    params.shieldingSolidAngle = cval("SHIELDING_SOLID_ANGLE");
  }

  void setParams() {
    params.mergeDistance = 0.03f;
    params.stepDistance = 0.02f;
    params.r_groundTissue = 1.75f;
    params.r_veins = 1.f;
    params.veinFullGrownAge = 250.f;
  }

  void loadReceptacle() {
    auto loaded = modelling::readMeshCageDeformationFromBinaryFile(
        "./bellis_contours_Mik_01_12_23.raw");
    //"./gerbera_contours.raw");
    if (loaded) {
      growingReceptacle.reload(std::move(*loaded));
      std::cerr << "Loaded contours: " << growingReceptacle.size() << '\n';
    }
  }
  void loadPhyllotaxis() {
    orderedPrims =
        phyl::readTimestampPrimordiaFromFile("./bellis_stable.prims");
    // phyl::readTimestampPrimordiaFromFile("./gerbera_initial.prims");

    // for (auto const &p : orderedPrims) {
    //  std::cerr << "PRE " << p << std::endl;
    //}

    orderedPrims = phyl::fixOutOfOrder(orderedPrims);
    // for (auto const &p : orderedPrims) {
    //  std::cerr << "POST " << p << std::endl;
    //}
  }

  //
  // SETUP
  //
  void init() override {
    frame = 0;
    randutil::seed(137);

    setParams();
    rereadPanel();

    loadReceptacle();
    loadPhyllotaxis();

    makeStemBase();
  }

  //
  // STEP SIMULATION
  //
  void step() override {

    rereadPanel();

    std::cerr << "Bellis new Step " << StepNo()
              << " time: " << growingReceptacle.currentTime() << '\n';

    // if (growingReceptacle.currentTime() > 0.65)
    //  return;

    // Features:
    // - modelled phyllotaxis
    // - lateral/extra veins for trifurcation/dense bracts
    addPrimordia();

    // Features:
    // - grow new veins (just tips single shot)
    // - minimal resistance path from auxim maxima to veins
    // -    parameters: current trajectory, step size, merge radius, ...
    //
    // - sectoriality (optional)
    // - randomness (optional)
    // - age dependency (optional, dependent on resistivity as function of
    // age)
    //
    if (growingReceptacle.currentTime() > 0.1f) {
      auto veinTips = gatherVeinTips();

      for (auto &tip : veinTips) {
        tip->data.age = 0;

        auto newTip = growVeinInitialOrthogonalToReceptical(
            tip, veins, params.stepDistance, params.veinJitterStrength,
            params.veinJitterFrequency);

        if (tip->data.type == MRRodkaew::NodeType::FLORET) {
          auto node = tip;
          while (node != newTip) {
            node->data.type = MRRodkaew::NodeType::FLORET;
            node->data.baseResistance = 100.f;
            node = node->parents()[0];
          }
          newTip->data.type = MRRodkaew::NodeType::ADAXIAL;
        }

        tip = newTip;
      }

      vasc::growMinimalResistancePathVeins(
          veins, veinTips,       //
          params.r_groundTissue, //
          params.mergeDistance,  //
          params.stepDistance, params.veinJitterFrequency,
          params.veinJitterStrength, params.sectorialPropensity,
          params.shieldingSolidAngle);
    }

    // Optional Features:
    // -reticulation
    if (params.useReticulation) {
      auto newNodes = vasc::processReticulateNodes(
          veins,                        //
          reticulateNodesToBeProcessed, //
          reticulateNodes,              //
          params.mergeDistance,         //
          params.stepDistance, params.veinJitterFrequency,
          params.veinJitterStrength, params.sectorialPropensity,
          params.reticulationPropensity, params.r_groundTissue,
          params.reticulationAgeTrigger);

      for (auto node : newNodes) {
        subReticulateNodes.push_back(node);

        if (node->data.type == MRRodkaew::NodeType::FLORET)
          node->data.type = MRRodkaew::NodeType::ADAXIAL;
      }
    }

    // -diffuse auxin and veins: space colonization
    // if (growingReceptacle.currentTime() > 0.48f)
    //  addDiffuseAuxin();
    // if (growingReceptacle.currentTime() > 0.52f) {
    //  auto veinTips = vasc::growDiffuseVeins(
    //      veins, diffuseAuxin, params.r_groundTissue, //
    //      params.mergeDistance,                       //
    //      params.stepDistance * 0.35f,                // TODO fix
    //      params.veinJitterFrequency, params.veinJitterStrength,
    //      params.sectorialPropensity, params.killDistance);
    //  for (auto tip : veinTips) {
    //    tip->data.type = MRRodkaew::NodeType::ADAXIAL;
    //    tip->data.baseResistance = 1.f;
    //  }
    //}

    // Features:
    // - Volume of Revolution
    // - Data-driven cage based deformation
    // - subdivide vein segments that get too long (merge as well, but less
    // likely)
    receptacleGrowthAndVeinDeformation();

    // Features:
    // -Effective resistance based on vein age, resistivity, and length
    auto t = growingReceptacle.currentTime();
    vasc::stepVascularAge(veins,                 //
                          params.r_groundTissue, //
                          params.r_veins,        //
                          params.veinFullGrownAge *
                              func(vein_age_fullgrown_f, t));

    int ageFullyGrown = params.veinFullGrownAge * func(vein_age_fullgrown_f, t);
    for (auto node : veins) {

      auto age = math::map(node->data.age, 0, ageFullyGrown, 0.f, 1.f);
      // node->data.resistanceAtNode = func(r_vein_node_f, age);

      node->data.resistanceAtNode =
          math::lerp(params.r_groundTissue, 1.f,
                     math::smoothstep(node->data.age, 0.f, ageFullyGrown));
    }

    vasc::acropetal_recalculateVascularResistances(veins, stemBase);

    // For visualization
    // put here for continuous panel modeling, else put in step
    recalculateVascularWidths();

    frame++;
  }

  //
  // DISPLAY SIMULATION
  //

  void drawBoundingBox(geometry::Curve2f controlPoints) {
    {
      auto maxBox = math::Vec3f{-9999.f, -9999.f, -9999.f};
      auto minBox = math::Vec3f{9999.f, 9999.f, 9999.f};

      auto cp = math::vec3f::from(controlPoints.front());
      auto c = cp;
      maxBox = math::maxComponents(c, maxBox);
      minBox = math::minComponents(c, minBox);
      c = cp;
      c.x *= -1;
      maxBox = math::maxComponents(c, maxBox);
      minBox = math::minComponents(c, minBox);
      c = cp;
      c.z = c.x;
      c.x = 0;
      maxBox = math::maxComponents(c, maxBox);
      minBox = math::minComponents(c, minBox);
      c = cp;
      c.z = -c.x;
      c.x = 0;
      maxBox = math::maxComponents(c, maxBox);
      minBox = math::minComponents(c, minBox);

      nproduce SetColor(1) SetWidth(0.01) MoveToV3(cp);

      for (int i = 1; i < controlPoints.size(); ++i) {
        cp = math::vec3f::from(controlPoints[i]);

        auto c = cp;
        maxBox = math::maxComponents(c, maxBox);
        minBox = math::minComponents(c, minBox);
        c = cp;
        c.x *= -1;
        maxBox = math::maxComponents(c, maxBox);
        minBox = math::minComponents(c, minBox);
        c = cp;
        c.z = c.x;
        c.x = 0;
        maxBox = math::maxComponents(c, maxBox);
        minBox = math::minComponents(c, minBox);
        c = cp;
        c.z = -c.x;
        c.x = 0;
        maxBox = math::maxComponents(c, maxBox);
        minBox = math::minComponents(c, minBox);

        // nproduce LineToV3(cp);
      }

      auto dim = maxBox - minBox;
      auto padDim = dim;
      dim += math::Vec3f{1.f, 1.f, 1.f};
      padDim -= (dim);
      minBox += 0.5 * padDim;

      // Draw
      // nproduce MoveToV3(minBox) SetHead(1, 0, 0, 0, 1, 0) F(dim.x) Up(90)
      //    F(dim.y) Up(90) F(dim.x) Right(90) F(dim.z) Right(90) F(dim.x)
      //    Up(90)
      //        F(dim.y) Up(90) F(dim.x);

      // invisible
      nproduce MoveToV3(minBox) SetHead(1, 0, 0, 0, 1, 0) f(dim.x) Up(90)
          f(dim.y) Up(90) f(dim.x) Right(90) f(dim.z) Right(90) f(dim.x) Up(90)
              f(dim.y) Up(90) f(dim.x);
    }
  }

  void drawAngleField() {
    /*  auto const &currentCage = growingReceptacle.currentCage();
      auto const &fields = currentCage.controlPointFields;
      auto const &mesh = currentCage.mesh;

      std::vector<float> angles(currentCage.controlPoints.size(), 0.f);
      for (int id = 0; id < angles.size(); ++id) {
        angles[id] = 360.f * (float(id) / (angles.size() - 1));
      }

      for (int id = 0; id < mesh.vertexCount(); ++id) {
        float angle{0.f};
        for (int cpID = 0; cpID < currentCage.controlPoints.size(); ++cpID) {
          angle += angles[cpID] * fields[cpID][id];
        }

        auto const &v = mesh.vertices()[id];
        auto dir = math::rotateAroundZ(math::axis3::X, angle);
        drawLine(v, v + math::vec2f::from(dir) * 0.1f);
      }*/
  }

  void drawCageControlPointWeights(int cpID) {
    // auto const &currentCage = growingReceptacle.currentCage();

    // auto w = currentCage.controlPointFields[0];
    // drawMesh(currentCage.mesh, // triangle mesh
    //         w, [](float w) -> int {
    //           return math::map(std::pow(w, 0.2f), 0.f, 1.f, 64, 255);
    //         });
  }

  void drawContour() {
    auto const &meshCage = growingReceptacle.currentCage();
    auto receptacle = geometry::subCurvePartition(
        meshCage.controlPoints, 0, meshCage.controlPoints.size() - 20);

    receptacle = geometry::reverse(receptacle);
    receptacle = geometry::subdivide(receptacle, 2);

    auto offset{-0.1f};
    geometry::Curve2f::point_store points;
    points.push_back(receptacle.front() + math::Vec2f(-offset, 0.f));
    for (uint32_t index = 1; index < receptacle.size() - 1; ++index) {
      auto a = receptacle[index - 1];
      auto b = receptacle[index];
      auto c = receptacle[index + 1];

      auto ab = b - a;
      auto bc = c - b;

      ab = perpendicularTo(ab);
      ab.normalize();

      bc = perpendicularTo(bc);
      bc.normalize();

      auto dir = normalized(ab + bc);
      points.push_back(b + offset * dir);
    }
    points.push_back(receptacle.back() + math::Vec2f(0.f, -offset));
    receptacle.points() = points;

    auto flip = receptacle;
    for (auto &p : flip) {
      p.x *= -1;
    }

    auto weight = func(view_size_f, growingReceptacle.currentTime());
    auto viewHeight = 1.f * weight + 4.f * (1.f - weight);

    nproduce SetWidth(0.05f) SetColor(10);
    drawCurve_GC(receptacle);
    drawCurve_GC(flip);

    /*
    receptacle = geometry::offsetContour(receptacle, 0.2f);
    receptacle = geometry::reverse(receptacle);

    nproduce SetColor(5) SetWidth(0.1);
    drawMirroredDashedCurve(receptacle, 0.2f);
    */
  }

  void drawClosestNeighbours() {
    nproduce SetWidth(0.1f);
    for (auto node : veins) {
      int nID = 0;
      if (node->data.closestNeighs.empty())
        continue;

      for (auto const neigh : node->data.closestNeighs) {
        nID++;
        nproduce SetColor(nID);

        if (viz.neighbourCombo == 1 && nID == 3)
          continue;
        if (viz.neighbourCombo == 2 && nID == 2)
          continue;
        if (viz.neighbourCombo == 3 && nID == 1)
          continue;

        drawLine(node->data.position, lerp(node->data.position, neigh->position,
                                           viz.neighbourLength));
      }
    }
  }

  void drawOldestNeighbours() {
    nproduce SetWidth(0.1f);
    for (auto node : veins) {
      int nID = 0;
      if (node->data.oldestNeighs.empty())
        continue;
      // drawSphere(node->data.position, 0.05f);
      for (auto const neigh : node->data.oldestNeighs) {
        // if (node->data.oldestNeighs.size() > 2) {
        //  auto neigh = node->data.oldestNeighs[1];
        nID++;
        nproduce SetColor(nID);

        if (viz.neighbourCombo == 1 && nID == 3)
          continue;
        if (viz.neighbourCombo == 2 && nID == 2)
          continue;
        if (viz.neighbourCombo == 3 && nID == 1)
          continue;

        drawLine(node->data.position, lerp(node->data.position, neigh->position,
                                           viz.neighbourLength));
      }
    }
  }

  void drawVasculature_Edges(
      std::function<int(Vasculature::node_ptr)> const &colorOf) {

    for (auto node : veins) {
      if (!node->hasParents())
        continue;
      auto cID = colorOf(node);
      auto w = veinWidth(node);
      if (node->isDisconnected())
        continue;
      nproduce SetWidth(w) SetColor(cID);
      drawSphere(node->data.position, w * 0.5f);
      for (auto child : node->children()) {
        auto cID = colorOf(child);
        auto w = veinWidth(child);

        nproduce SetWidth(w) SetColor(cID);
        drawLine(node->data.position, child->data.position);
      }
      if (node->isLeaf()) {
        auto w = veinWidth(node);
        auto cID = colorOf(node);
        nproduce SetColor(cID);
        drawSphere(node->data.position, w * 0.5f);
      }
    }
  }

  void drawVasculature_GeneralizedCylinders(
      std::function<int(Vasculature::node_ptr)> const &colorOf) {

    for (auto stem : stemBase) {
      if (!stem->hasChildren())
        continue;

      std::queue<Vasculature::node_ptr> branches;

      for (auto child : stem->children()) {
        branches.push(child);
      }

      while (!branches.empty()) {
        auto node = branches.front();
        branches.pop();

        auto next = node;
        node = node->parents()[0];
        auto dir = math::normalized(next->data.position - node->data.position);
        auto head = V3f(dir.x, dir.y, dir.z);

        auto cID = colorOf(next);
        auto w = veinWidth(next);

        auto type = next->data.type;

        nproduce                          //
        SB()                              //
            MoveToV3(next->data.position) //
            SetWidth(w)                   //
            SetColor(cID)                 //
            SetHead3f(head)               //
            StartGC();                    //

        node = next;
        while (node->hasChildren()) {

          int childID = 0;
          next = node->children()[childID++];
          // while (next->data.type != type && childID <
          // node->children().size()) {
          //  branches.push(next);
          //  next = node->children()[childID++];
          //}

          for (int i = childID; i < node->children().size(); ++i) {
            branches.push(node->children()[i]);
          }
          // if (next->data.type != type) {
          //  break;
          //}

          auto dir =
              math::normalized(next->data.position - node->data.position);
          auto head = V3f(dir.x, dir.y, dir.z);

          auto cID = colorOf(next);
          auto w = veinWidth(next);

          nproduce SetColor(cID) SetWidth(w) SetHead3f(head)
              LineToV3(next->data.position) PointGC();

          node = next;
        }

        nproduce EndGC() EB();

        if (node->isLeaf()) {
          auto w = veinWidth(node);
          auto cID = colorOf(node);
          nproduce SetColor(cID);
          drawSphere(node->data.position, w * 0.5f);
        }
      }
    }
  }

  void draw() override {

    drawBoundingBox(growingReceptacle.currentCage().controlPoints);

    if (viz.render.dashedContour)
      drawContour();

    if (viz.render.closestNeighbours)
      drawClosestNeighbours();

    if (viz.render.oldestNeighbours)
      drawOldestNeighbours();

    std::function<int(Vasculature::node_ptr)> colorOf;
    colorOf = [&](Vasculature::node_ptr node) { return colorID(node); };
    if (viz.render.style == Viz::Render::Style::AgeResistance) {
      colorOf = [=](Vasculature::node_ptr node) {
        auto r = math::clamp(node->data.resistanceAtNode, params.r_veins,
                             params.r_groundTissue);
        return math::map(r,                                     //
                         params.r_groundTissue, params.r_veins, //
                         64, 255);
      };
    }
    if (viz.render.style == Viz::Render::Style::EffectiveResistance) {
      auto maxVal = 0.f;
      for (auto node : veins) {
        if (!node->hasChildren())
          continue;

        if (node->data.effectiveResistance > maxVal)
          maxVal = node->data.effectiveResistance;
      }
      colorOf = [=](Vasculature::node_ptr node) {
        return math::map(node->data.effectiveResistance, 0.f, maxVal, 64, 255);
      };
    }
    // drawVasculature_GeneralizedCylinders(colorOf);
    drawVasculature_Edges(colorOf);

    // draw primordia
    // nproduce SetColor(viz.color.primordia);
    // for (auto prim : primordia) {
    //  auto w = viz.size.primordia *
    //           (2.f + math::smoothstep(float(prim->data.age), 0.f, 500.f)) *
    //           prim->data.width * 0.5;
    //  drawSphere(prim->data.position, w);
    //}

    // draw florets
    // for (auto floret : florets) {
    //  auto w = viz.size.floret * veinWidth(floret) * 0.5f;

    //  auto colorID = viz.color.floretPrim;
    //  if (floret->data.position.z < -viz.backsideDepth)
    //    colorID = viz.color.backside;

    //  nproduce SetColor(colorID);
    //  drawSphere(floret->data.position, w);
    //}

    //// draw bracts
    // for (auto bract : bracts) {
    //  auto w = viz.size.bract * veinWidth(bract) * 0.5f;

    //  auto colorID = viz.color.bractPrim;
    //  if (bract->data.position.z < -viz.backsideDepth)
    //    colorID = viz.color.backside;

    //  nproduce SetColor(colorID);
    //  drawSphere(bract->data.position, w);
    //}

    //// draw stem base
    // nproduce SetColor(viz.color.stem);
    // for (auto sink : stemBase) {
    //  auto w = viz.size.stem * sink->data.width * 0.5f;
    //  if (sink->data.effectiveResistance > 0.f)
    //    nproduce SetColor(255);
    //  else
    //    nproduce SetColor(6);
    //  drawSphere(sink->data.position, w);
    //}

    // draw reticulation nodes
    // nproduce SetColor(viz.color.stem);
    // for (auto node : reticulateNodes) {
    //  // auto w = viz.size.stem * node->data.width * 0.5f;
    //  auto w = 0.01f;
    //  drawSphere(node->data.position, w);
    //}

    // draw reticulation nodes to be processed
    // nproduce SetColor(viz.color.floret);
    // for (auto node : reticulateNodesToBeProcessed) {
    //  // auto w = viz.size.stem * node->data.width * 0.5f;
    //  auto w = 0.01f;
    //  drawSphere(node->data.position, w);
    //}

    nproduce SetWidth(0.f);
    drawCageControlPointWeights(50);

    nproduce SetWidth(0.1f);
    drawAngleField();

    // draw diffuse auxin
    if (viz.render.diffuseAuxin) {
      nproduce SetColor(viz.color.auxin);
      for (auto const &auxinPoint : diffuseAuxin) {
        drawSphere(auxinPoint, 0.1f);
      }
    }

    // dbg draw just nodes
    // for (auto node : veins)
    //{
    //    auto cID = colorOf(node);
    //    auto w = veinWidth(node);
    //    nproduce SetColor(cID);
    //    drawSphere(node->data.position, w);
    //}
  }
};

//
// LPFG Main loop
//
// clang-format off
derivation length: 1;

module Draw();

Axiom: Draw();

Start: {
	// should be std::make_unique however preproc.sh only works with c++11 by default
	// If you change preproc.sh to use c++17 everything works fine with std::*** 17 functions
	model = util::make_unique<BellisVasculature>();

    model->init();
}

StartEach: {
    model->step();
}

interpretation:

Draw(): {
	model->draw();
	produce;
}
// clang-format on
