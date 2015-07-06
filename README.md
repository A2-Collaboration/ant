ant
===

Just another analysis toolkit `ant`, based on GoAT tree output from
AcquRoot.


## GoAT File Reader

Almost implemented.

### Notes
 * Output module missing
 * Physics manager missing

#### Tree Manager
TODOs:
 * remove tree is module that requested it fails to init. Right now only the module is removed, leaving the tree active
 * keep track of associated branches -> Allow modules to share trees
 * unset branch addresses on module destruction? -> Possible (probable) segfault here, if module terminates after setting a branch and tree wants to read into it
 * unset unused branches? -> speed?

#### Particles
Currently negelecting energy and angle inforamtion from particle trees
and just using the referenced track. This info is a dublicate of the
track data anyway. Ohterwise no idea what to do with it.

#### Detector Hits
Ant data structure missing.

#### Pluto
Works, but no decay tree inforamtion. PParticle::GetDauther() always returns nullptr.
