# Changelog

&nbsp;

A complete list of all changes can be found here:

https://github.com/A2-Collaboration/ant/commits/master


## Noteworthy Changes

Below is an incomplete list of changes which might be of general interest.

&nbsp;

 * Move `Matches()` to base Setup and specify the time range instead via `SetTimeRange(start, end)`; start and end date can now be queried
 * Add support for 1D histograms with a variable bin width to `HistogramFactory` (see also `VarBinSettings` and `VarAxisSettings`)
 * ...


## v1.2

 * Take increased Mk2 buffer size into account, implemented new method to automatically determine it
 * Fix CMake to work with new official pluto6 fork
 * Try to avoid Clusters with `NaN` times if they're marked as `BadTDC`, use time from neighbouring crystal of the cluster instead
 * Helper functions added to switch off Tagger range and get a list of BaF2 and PbWO4 channels
 * More Setups added
 * Minor improvements and changes

 * Ant-mcgun: Add option to generate particles with a flat theta distribution

 * Tests: ROOT6 integrated in CI with Travis


## v1.1

Finalized ROOT6 migration, minor improvements
 * Fix some warnings related to the ROOT6 migration
 * Fix linker problem on some recent OS systems
 * Fix warnings when compiling with Clang 6.0

 * Add some easy methods for manipulating and modifying hstack objects, like Rebin, Add, Multiply, and Divide
 * Code moved from A2-Collaboration-dev to A2-Collaboration

 * Update documentation


&nbsp;

&nbsp;

All changes up to the tagged release 1.0 are omitted. Starting with the release of version 1.0 Ant could be used with ROOT6 and the above changelog reflects notable changes to the codebase.
