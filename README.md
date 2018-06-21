[![Build status master branch](https://travis-ci.org/A2-Collaboration/ant.svg?branch=master)](https://travis-ci.org/A2-Collaboration/ant)
[![Codecov status master branch](https://codecov.io/github/A2-Collaboration/ant/coverage.svg?branch=master)](https://codecov.io/github/A2-Collaboration/ant?branch=master)
[![Coverity Scan build status](https://scan.coverity.com/projects/7274/badge.svg)](https://scan.coverity.com/projects/a2-collaboration-dev-ant)
[![Gitter](https://img.shields.io/gitter/room/A2-Collaboration-dev/Lobby.svg)](https://gitter.im/A2-Collaboration-dev/Lobby)

ant
===

![A2 ant logo](doc/a2-ant-logo.png?raw=true)

Just another **AN**alysis **T**oolkit `ant`, which can read from many
input sources with minimal user intervention and let's you create
Physics analysis within minutes.

Please see also the automatically generated
[Doxygen pages](http://a2-collaboration-dev.github.io/ant/).


## Dependencies
  * C++11 (gcc >4.9.2 should work, but >6.0 recommended, clang is also supported)
  * cmake >3.0
  * libgsl-dev, liblzma-dev (before installing root)
  * [CERN ROOT5](https://root.cern.ch/) (ROOT6 is supported as well)
  * [Pluto](https://github.com/A2-Collaboration/gsi-pluto) (if you use ROOT6, you have to use the pluto6 branch)
  * [APLCONpp](https://github.com/A2-Collaboration/APLCONpp)
  * doxygen (optional)

## Installation

**Note:** A comprehensive installation guide can be found in the [Wiki](https://github.com/A2-Collaboration/ant/wiki/How-to-Install).

Please make sure that you fulfill the dependencies. Install
[ROOT](https://root.cern.ch/building-root) and Pluto preferably in
`/opt/root` and `/opt/pluto` to ensure auto-detection by ant's CMake.

Next, you need to get the
[APLCON C++ wrapper](https://github.com/A2-Collaboration/APLCONpp).
The easiest way is to clone the repository relative to your `ant`
directory at `../APLCONpp`, then create the build directory at
`../APLCONpp/build`. This way CMake will automatically detect it.

See also the detailed steps described in the corresponding [Wiki section](https://github.com/A2-Collaboration/ant/wiki/How-to-Install#install-aplcon-with-c-wrapper).

Now you should be able to compile the Ant framework.
Therefore clone this respository, either directly from
https://github.com/A2-Collaboration/ant.git or you may want to fork it.
As mentioned above, the installation works best when the APLCONpp and ant folder reside in the same directory.
Inside your ant directory create a build direcory:

`mkdir build && cd build && cmake ..`

Then start the parallel compilation for example on a QuadCore machine with `make -j5`.
You may want to add your `ant/build/bin` directory to your `$PATH` variable.

## Parallel Processing

Ant is designed to run single threaded. This avoids a lot of programming
problems and running it on the computation clusters is simpler and more
efficient, since you only need to allocate one thread per job. For processing
data we recommend to run an Ant process on each input file and then merge the
results afterwards using `Ant-hadd`, `Ant-chain`, or ROOTs `hadd` tool (but prefer
`Ant-hadd --native` for custom data type support).

There is also no builtin option to run over multiple input files in
one go. This should be handled by external tools like GNU `parallel`, or
`AntSubmit` on a cluster (see also `--no_qsub` option), or your shell.
See below for an `AntSubmit` quick start guide.

## Quick start guides

Check the Wiki to learn about the basic usage of [Ant](https://github.com/A2-Collaboration/ant/wiki/Running-Ant) itself
or how to run Ant tools on [blaster](https://github.com/A2-Collaboration/ant/wiki/Blaster) using the provided job submission scripts.

### Ant MC tools
Ant comes with a few tools to generate MC data.
It can generate photoproduction and decays with `Ant-pluto`, a frontend utilizing Pluto for A2 physics (includes the Tagger),
shoot particles randomly in all directions using `Ant-mcgun`,
or simulate a complete cocktail of various photoproduction and decay channels according to their cross sections with `Ant-cocktail`.

#### Example: Pluto Decay
To use Pluto to simulate, for example, the omega ---> pi0 gamma do:
```
Ant-pluto --reaction "p omega [ pi0 g ]" --Emin 1400 --Emax 1600 --numEvents 10000 -v 2 -o sim.root
```
This will generate 10k events in the incident photon energy range 1400 MeV to 1600 MeV, saving unstable particles.
The pi0 will decay into different channels according to the Pluto database.
Use the option `--no-unstable` if you don't want to store intermediate particles and `--no-bulk` to disable bulk decays using the Pluto database.

#### Example: Random Gun
Shoot 1000 protons into TAPS, 1 proton/Event, 0 to 1 GeV
```
Ant-mcgun --numEvents 1000 --particle p --theta-max 25 --Emax 1000
```

## Troubleshooting

  * In case CMake is not able to locate your Pluto installation, you can provide the environment variable `$PLUTOSYS` to tell CMake where to find it. If you installed Pluto inside your home directory, `~/pluto` or `~/src/pluto`, or placed it in `/opt/pluto`, the make process may have failed. Please make sure you ran `make` in your Pluto directory with a proper ROOT installation.
  * If you're using gcc version 5.x and experiencing build errors within ROOT generated dictionary files (i.e. redeclaration of `struct __xfer_bufptrs`), please update to a more recent ROOT version. As of November 2015, you need to clone the git branch which includes the patches. To do so:
    * get the source with `git clone -b v5-34-00-patches https://github.com/root-mirror/root.git`
    * `cd` in the cloned directory
    * create a new folder for the build process, e.g. `mkdir build_dir`
    * `cd` in it and run `cmake .. && make -jN`, replace `N` with the number of threads which should be used
    * set your `$ROOTSYS` accordingly
  * gdb version 7.7 crashes in combination with cereal, so use version 7.10
  * If you are using ROOT6 and try to store own classes, even derived from ROOT classes, nested in a STL container in a tree, like `std::vector<ant::TSimpleParticle>`, you might get weird errors starting with `Fatal in <TBranchElement::InitializeOffsets>` by Ant. This is a known bug introduced in ROOT6 and reported here: [ROOT-9450](https://sft.its.cern.ch/jira/browse/ROOT-9450). If you get those errors, try to avoid nesting the objects in a container and store them individually instead as long as the bug is not fixed.
  * If using the newest GCC release, like 7.x as of June 2017 or higher in the future (usually ArchLinux users), use the patched ROOTv5 version available [here](https://github.com/A2-Collaboration/root/commits/v5-34-00-patches-A2), which contains patches based on the bug report [ROOT-8180](https://sft.its.cern.ch/jira/browse/ROOT-8180). In the meantime, those patches have been merged upstream.
  * In case you encounter weird test errors (garbage output of corrupted TFiles `TKey::ReadObjWithBuffer msg='Unknown class foDoubly linked listZL`), make sure you recompiled Pluto properly with the identical `$ROOTSYS` as Ant (in particular pay attention to ROOT debug/release builds).
  * Starting with clang 5.0, there's a problem with compiling the third party tool cereal. The problem is [known and reported](https://github.com/USCiLab/cereal/issues/439). Please use a clang version < 5.0 as long as this issue with RapidJSON isn't fixed.

# Contributing

Please read the following sections if you want to contribute to this project.
Always make sure to cover your code with tests, see the [wiki](https://github.com/A2-Collaboration/ant/wiki/Tests).
Run `make build_and_test` before pushing.

## Coding Style
  * Indentation: 4 spaces, no tabs anywhere
  * Ant-Codesytle defined in [doc/Ant-Codestyle.xml](doc/Ant-Codestyle.xml) (import to QtCreator)
  * `#include` statements:
    * header file for this .cc file first
    * then grouped by ant, ROOT, STL, others
    * one line between groups
    * each group ordered alphabetically

## Library organization

The calibration modules specify physics classes, usually below
`src/analysis/physics/calibration`, which produce the histograms used
in the GUI part of the module. The analysis classes can access all
information organized in `expconfig`. The `slowcontrol` libraries are
divided to ensure proper initialization order of static registries. 

<img src="doc/library-dependencies.png">

The solid arrows mean "links to", whereas the dashed arrows means
"includes only".

Or, automatically generated with
```bash
cd build
cp ../doc/CMakeGraphVizOptions.cmake .
cmake --graphviz=CmakeGraphViz ..
dot -Tpng CmakeGraphViz.Ant -o ../doc/library-dependencies-autogenerated.png
rm CmakeGraphViz*
```
<img src="doc/library-dependencies-autogenerated.png">

## Development

The following items are still to-do:

  * Implement the wire chamber detector (hard), or the conventional
    tagger ladder including magnetic field energy calibration
    (easier)
  * ~~Implement Mk1 unpacker (many things already provided)~~ Done including GZ decompression. Scaler buffer decoding still WIP.
  * ~~Implement EPICS reader, and some more slow control variables~~ Done for tagging efficiency at least.  

### Detector Type Mapping to Goat/Acqu

| Ant       | Goat |
|-----------|------|
| CB        | NaI  |
| PID       | PID  |
| MWPC0     | MWPC |
| MWPC1     |      |
| TAPS      | BaF2 |
|           | PbWO4|
| TAPSVeto  | Veto |

All other Ant detectors to not have a representation in Goat/Acqu.

### Data flow

The physics classes analyse `TEvent`s provided by different sources or
amenders. The main source of events is the AntReader, which itself is
either fed by some unpacker or by already unpacked and possibly
reconstructed `treeEvents`. Additionally, the SlowControlManager needs
to do proper buffering to make physics classes easy to implement.

<img src="doc/dataflow.png">

## External Components

Have a look at those very nice projects, which are used here:

  * [Easylogging++](http://easylogging.muflihun.com/)
  * [Catch](https://github.com/philsquared/Catch) framework for unit-tests, test-driven development. See [the test/ subdirectory](test/).
  * [TCLAP - Templatized C++ Command Line Parser](http://tclap.sourceforge.net)
  * [cereal](http://uscilab.github.io/cereal/) for [ant::TEvent](src/tree/TEvent.h) serialization into ROOT TTree
