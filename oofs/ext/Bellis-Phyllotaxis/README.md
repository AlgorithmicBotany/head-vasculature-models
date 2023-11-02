# A model of phyllotaxis for Bellis (common daisy)

![Bellis sequence](bellis_sequence_phyllotaxis.png)

This model is a slight modification of the gerbera model from Zhang et al. (2021).
More precisely, it is a modified version of Model 3 from that paper, which can be downloaded from here:
http://algorithmicbotany.org/papers/gerbera2021.html

The key differences are:

- The curves defining the shape of the head, which were changed to match longitudinal sections of Bellis, and

- The function defining the threshold distance for primordium insertion, which was changed to model the difference in sizes
between Bellis and gerbera primordia.

The vlab object contains three seperate lpfg models:

- The phyllotaxis model, which simulates patterning on the receptacle. Its components are listed under the sub-heading
`PHYLLOTAXIS MODEL` in the vlab object's menu.

- The contour model, which is used to define the profile curves and landmarks of the receptacle growth model. Its
components are listed under the sub-heading `LONGITUDINAL SECTIONS`.

- The receptacle growth model, which simulates the growth of the receptacle by interpolating the profile curves from the contour model
and defines the knots of the B-Spline basis functions used in the interpolation. Its components are listed under the sub-heading
`RECEPTACLE MODEL`.

The user-defined functions and contours that are common to all three models are accessible from the `GALLERY` sub-heading
in the vlab object. If the landmarks are modified in the `contours_markers.cset` file, then you must regenerate the
`markers.h` file by running `generate` under the `MARKERS` sub-heading.

## Representation of the head

According to the current understanding of phyllotaxis at the molecular level, the key patterning processes 
take place in the epidermis of the meristem.  Consequently, we only represent
the (curved) surface of the plant.  The surface is parametrized by:

- The azimuth angle `phi`, measured (in the anti-clockwise direction) from the direction of the x axis, and 

- The curvilinear coordinate `s`, defined as the arc-length distance of points on the surface, measured from the base of the stem.  

Some aspects of the simulation, in particular calculation of distances between primordia and visualization, are conveniently 
performed using cylindrical and Cartesian coordinates.  Function `SyncCoordinates` converts the `phi`-`s` coordinates into
cylindrical and Cartesian coordinates

## Modeling growth

Plant growth is modeled descriptively, by explicitly specifying the changes in shape and size of the receptacle.
 
- Changes in the shape of the receptacle are modeled by interpolating between a number of curves that represent 
its profile and different stages of development. Both curves are defined graphically using the lpfg contour editor (cuspy)
and accessed from lpfg using its built-in function `curveXYZ`, which returns the xyz coordinates of a point on the curve given 
its arc-length coordinate. The interpolation is performed by computing a weighted average 
of the points with the same value of parameter `ns`. The weight is a variable `Progress`, which ranges from 0 to 1 and
represents the progress of growth over time.  

- Superimposed on the changes in shape are changes in size, modeled by uniformly scaling the profile curve by factor 
`CURVE_SCALE`.

## Placement of primordia

The model places primordia on a growing meristem/receptacle accrding to the Hofmeister/Snow and Snow rule interpreted as:
"place a new primordium on the active ring when and where there is enough space for it."
Primordia are represented as spheres with a given center and radius.  The distance between primordia is measured
as the sum of their radii.  With this interpretation, an incipient primordium is placed on the active ring 
when its distance from all previously placed primordia is greater than a threshold value. 
The distances are measured as the difference in angular position and in arc-length position.   

Primordia boundaries are represented as a set of points on the surface of the receptacle, with a given distance `rho0`
from the primordium center.  As the receptacle grows, both the
size and shape of the primordia embedded in it may thus change.  

A primordium is described in three different coordinate systems:
- arc-lengh-azimuth (used to simulate propagation of primordia on the surface of the meristem)
- cylindrical (intermediate step, before conversion to Cartesian)
- Cartesian (used to measure distances, and display results)


## Reference
```
@article{2021ZhangPNAS,
       author = {{Zhang}, Teng and {Cieslak}, Mikolaj and {Owens}, Andrew and {Wang}, Feng and {Broholm}, Suvi K. and {Teeri}, Teemu H. and {Elomaa}, Paula and {Prusinkiewicz}, Przemyslaw},
        title = "{Phyllotactic patterning of gerbera flower heads}",
      journal = {Proceedings of the National Academy of Science},
         year = 2021,
        month = mar,
       volume = {118},
       number = {13},
          eid = {e2016304118},
        pages = {e2016304118},
          doi = {10.1073/pnas.2016304118},
}
```




