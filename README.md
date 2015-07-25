# hyperelastic
Javascript simulations of hyperelastic materials submitted to growth.

#Adding a new geometry

These are the steps that were required to add the block geometry and the Block Border preset (and may be required to add another geometry)

to geometry.js
* add makeBlock()
* add blockTetra()
* add blockVind()

to growth.js
* add growBlockBorderInstantaneous()

to presets.js
* add a preset demonstrating the use of the geometry (optional)

to index.html
* add preset to the list of presets (if a preset had been added)
* add to the switch of handled url arguments  (if a preset had been added)
* add to the getPresetParams function (required, even if no new preset has been added)

to simulation.js
* add geometry to switch in initSimulation()
* add growth function to switch in initSimulation()

to display.js
* add block_deformationColor()

nothing to add to mechanics.js or to algebra.js