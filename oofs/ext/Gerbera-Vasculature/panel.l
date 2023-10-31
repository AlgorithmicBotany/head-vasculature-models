target: Lsystem
panel name: L-system
background: 0,0,0
size: 747 923
font: Arial,14,-1,5,50,0,0,0,0,0

type: PAGE
name: Configuration
color: 255,255,255
origin: 178 19
message: d NOTHING %d 1

type: SLIDER
name: Vein age full resistivity
colors: 255,255,255 255,0,0
origin: 53 318
min/max: 1 500
value: 13
message: d VEIN_AGE_FULLGROWN %d 1

type: LABEL
name: Visualization
color: 255,255,255
origin: 533 15

type: SLIDER
name: Vein scale
colors: 255,255,255 255,0,0
origin: 337 466
min/max: 0 10
value: 8
message: d VEIN_SCALE %d 100

type: SLIDER
name: Vein exponent
colors: 255,255,255 255,0,0
origin: 337 552
min/max: 1 20
value: 2
message: d VEIN_EXPONENT %d 10

type: SLIDER
name: Floret primordia scale /10
colors: 255,255,255 255,0,0
origin: 535 350
min/max: 0 40
value: 0
message: d FLORET_PRIM_SIZE %d 10

type: SLIDER
name: Bract Primordia Size /10
colors: 255,255,255 255,0,0
origin: 536 422
min/max: 0 40
value: 0
message: d BRACT_PRIM_SIZE %d 10

type: SLIDER
name: Primordia scale /10
colors: 255,255,255 255,0,0
origin: 535 282
min/max: 0 40
value: 0
message: d PRIM_SIZE %d 10

type: SLIDER
name: Vein tip width /10
colors: 255,255,255 255,0,0
origin: 337 638
min/max: 1 20
value: 20
message: d VEIN_TIP_WIDTH %d 10

type: SLIDER
name: Stem node scale /10
colors: 255,255,255 255,0,0
origin: 535 494
min/max: 0 20
value: 0
message: d STEM_NODE_SIZE %d 10

type: SLIDER
name: Render neighbours
colors: 255,255,255 255,0,0
origin: 74 440
min/max: 0 2
value: 0
message: d RENDER_NEIGHBOURS %d 1

type: LABEL
name: 0-None,1-Closest,2-Oldest
color: 255,255,255
origin: 60 494

type: SLIDER
name: Vein jitter frequency(%)
colors: 255,255,255 255,0,0
origin: 274 60
min/max: 0 100
value: 66
message: d VEIN_JITTER_FREQ %d 100

type: SLIDER
name: Vein jitter strength / 100
colors: 255,255,255 255,0,0
origin: 274 128
min/max: 0 100
value: 27
message: d VEIN_JITTER_STRENGTH %d 100

type: SLIDER
name: Sectorial propensity (%)
colors: 255,255,255 255,0,0
origin: 52 232
min/max: 0 100
value: 66
message: d SECTORIAL_PROPENSITY %d 100

type: SLIDER
name: Floret kick strength
colors: 255,255,255 255,0,0
origin: 53 60
min/max: 0 100
value: 38
message: d FLORET_KICK %d 100

type: SLIDER
name: Use reticulation
colors: 255,255,255 255,0,0
origin: 274 196
min/max: 0 1
value: 0
message: d USE_RETICULATION %d 1

type: SLIDER
name: Bract kick strength
colors: 255,255,255 255,0,0
origin: 52 146
min/max: 0 100
value: 16
message: d BRACT_KICK %d 100

type: SLIDER
name: Reticulation propensity (%)
colors: 255,255,255 255,0,0
origin: 273 264
min/max: 0 100
value: 100
message: d RETICULATION_PROPENSITY %d 100

type: SLIDER
name: Vein age fade in
colors: 255,255,255 255,0,0
origin: 526 575
min/max: 1 500
value: 5
message: d VEIN_AGE_FADE_IN %d 1

type: SLIDER
name: Adaxial vein scale /10
colors: 255,255,255 255,0,0
origin: 505 127
min/max: 0 40
value: 12
message: d ADAXIAL_VEIN_SCALE %d 10

type: SLIDER
name: Abaxial vein scale /10
colors: 255,255,255 255,0,0
origin: 506 199
min/max: 0 40
value: 16
message: d ABAXIAL_VEIN_SCALE %d 10

type: SLIDER
name: Floret vein scale /10
colors: 255,255,255 255,0,0
origin: 505 59
min/max: 0 40
value: 10
message: d FLORET_VEIN_SCALE %d 10

type: SLIDER
name: Render style
colors: 255,255,255 255,0,0
origin: 342 723
min/max: 0 2
value: 0
message: d RENDER_STYLE %d 1

type: SLIDER
name: Reticulation age trigger
colors: 255,255,255 255,0,0
origin: 266 340
min/max: 0 500
value: 254
message: d RETICULATION_AGE_TRIGGER %d 1

type: SLIDER
name: Render dashed
colors: 255,255,255 255,0,0
origin: 531 723
min/max: 0 1
value: 0
message: d RENDER_DASHED %d 1

type: SLIDER
name: Continuously reread panel
colors: 255,255,255 255,0,0
origin: 57 755
min/max: 0 1
value: 1
message: d CONTINUOUS_REREAD %d 1

type: SLIDER
name: Neighbour lengths
colors: 255,255,255 255,0,0
origin: 72 537
min/max: 0 100
value: 100
message: d NEIGHBOUR_LENGTH %d 100

type: SLIDER
name: Render neighbours combo
colors: 255,255,255 255,0,0
origin: 69 618
min/max: 0 3
value: 3
message: d NEIGHBOUR_COMBO %d 1

type: LABEL
name: 0-All,1-RG,2-RB,3-GB
color: 255,255,255
origin: 67 669

type: SLIDER
name: Render Diffuse Auxin
colors: 255,255,255 255,0,0
origin: 523 658
min/max: 0 1
value: 0
message: d RENDER_DIFFUSE_AUXIN %d 1

type: SLIDER
name: Sheilding (Pi)
colors: 255,255,255 255,0,0
origin: 339 816
min/max: 0 200
value: 5
message: d SHIELDING_SOLID_ANGLE %d 100

type: SLIDER
name: Backside cuttoff depth /100
colors: 255,255,255 255,0,0
origin: 63 842
min/max: 0 1000
value: 1000
message: d BACK_SIDE_AT_DEPTH %d 100

type: LABEL
name: Gerbera
color: 255,255,255
origin: 16 14

type: PAGE
name: Visualization
color: 255,255,255
origin: 163 5
message: p <pg#>
