# Changelog

&nbsp;

https://github.com/A2-Collaboration/ant/commits/master

&nbsp;

&nbsp;

&nbsp;

&nbsp;


Just kidding, below is a not complete list of changes which might be of general interest.

&nbsp;


 * Try to avoid Clusters with `NaN` times if they're marked as `BadTDC`, use time from neighbouring crystal of the cluster instead
 * ROOT6 integrated in continuous integration with Travis
 * ...


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
