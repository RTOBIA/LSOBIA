LSOBIA
======

__Large Scale Object Based Image Analysis__ (LSOBIA) is a module for
the [Orfeo Toolbox](https://www.orfeo-toolbox.org/) (OTB). It provides
several tools for object based, large scale remote sensing image,
analysis.

The module contains 5 OTBApplications:

* __LSSegmentation__ (Large Scale Segmentation) provides several
  methods to perform segmentation of very high resolution images.
* __LSSmallRegionsMerging__ (Large Scale Image Small Regions Merging)
  provides a method to perform small regions merging of very high
  resolution images.
* __Polygonize__ provides several methods to perform polygonization of
  high resolution images.
* __LSPolygonize__ (Large Scale Polygonize) is a distributed version
  of Polygonize. _work in progress_
* __ComputeAttributes__ computes several attributes of a vector file.

# Getting Started

## Dependencies

* [OTB](https://www.orfeo-toolbox.org/)
* The C++ library for the [Message Passing
  Interface](https://www.mpi-forum.org/) (MPI) (developed and
  validated with mpich-3.2)

Since LSOBIA is a module for the OTB, one need to install the OTB in
order to build it.

LSOBIA allows the distribution of the calculation on a computation
cluster, therefor it uses MPI for the communications between the
clusters.

## Building
LSOBIA can be built like any other [otb remote
module](https://wiki.orfeo-toolbox.org/index.php/How_to_write_a_remote_module)
You can build it either from within OTB's sources or outside it.

Don't forget to activate C++14 by setting the cmake parameter
_"CMAKE\_CXX\_FLAGS"_ to _"-std=c++14"_.

## Usage

LSOBIA can run on a single processor, or be distributed over multiple processors.

### Mono-processor execution

To run an application on a silgle processor, one need to call the
application using the OTB application launcher process:

```bash
otbApplicationLauncherCommandLine ${AppName} ${AppDirectory} ${Parameters}
```

* ${AppName} is the name of the application (ie: LSSegmentation or
  Polygonize).
* ${AppDirectory} is the path to the directory containing the compiled
  applications.
* ${Parameters} are the parameters for the application.

Exemple:

```bash
otbApplicationLauncherCommandLine LSSegmentation /home/me/bin/lsobia/lib/otb/applications -io.im inputimage.tif -io.out.dir /home/me/out -io.temp /tmp -algorithm baatz -algorithm.baatz.maxiter 45 -processing.memory 10000 -processing.nbproc 8 -processing.nbtilesperproc 2 -processing.writeimages on -processing.writegraphs off
```

### Multi-processor execution

Simply add "*mpirun -np ${NumberProcessor}*" before the previous
command. ${NumberProcessor} is the number of processors to be used for
the computation.

Exemple:

```bash
mpirun -np 8 otbApplicationLauncherCommandLine LSSegmentation /home/me/bin/lsobia/lib/otb/applications -io.im inputimage.tif -io.out.dir /home/me/out -io.temp /tmp -algorithm baatz -algorithm.baatz.maxiter 45 -processing.memory 10000 -processing.nbproc 8 -processing.nbtilesperproc 2 -processing.writeimages on -processing.writegraphs off
```

# Output Samples

## Baatz segmentation

We applied the LSSegmentation application with the Baatz algorith on
this image, and obtained the following label image as output. The
process used a single processor for the computation.

![input file](https://cloud.githubusercontent.com/assets/26165185/24074238/073c72b4-0c05-11e7-82e3-497a28d10db1.png) : 
![baatz-segmentation](https://cloud.githubusercontent.com/assets/26165185/24074247/3a061d62-0c05-11e7-9266-31a16d203afc.jpg)

# License
This project is licensed under the Apache License 2.0. Please see the
[LICENSE file](LICENSE) for legal issues on the use of the software.
