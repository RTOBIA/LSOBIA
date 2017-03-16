# LSOBIA

OTB Module providing Large Scale Object Based Image Analysis functionalities.

How to build it
==============
LSOBIA can be built like any other [otb remote module](https://wiki.orfeo-toolbox.org/index.php/How_to_write_a_remote_module)
You can build it either from within OTB's sources or outside it.

How to use it
============
LSOBIA provides an OTBApplication, LSSegmentation (Large Scale Segmentation).

```bash
  otbcli_LSSegmentation -io.im ${INPUT_IMAGE_PATH} -io.out ${OUTPUT_DIRECTORY} -algorithm baatz -algorithm.baatz.numitfirstpartial 1 -algorithm.baatz.numitpartial 1 -algorithm.baatz.stopping 10 -algorithm.baatz.spectralweight 0.05 -algorithm.baatz.geomweight 0.95 -processing.memory 2000 -processing.maxtilesizex 1000 -processing.maxtilesizey 1000 -io.temp ${TEMP_DIRECTORY} -processing.writeimages "on" -processing.writegraphs "on" -processing.aggregategraphs "on"
```

This produces an image output containing labels of the Baatz Segmentation algorithm.

Dependencies
===========
[OTB](https://www.orfeo-toolbox.org/), MPI

Licence
=======
Please see the license for legal issues on the use of the software.
