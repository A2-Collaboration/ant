ant-database
============

This folder is left blank and the corresponding .gitignore ensures
this. Every user should create symlinks where he wants the database
files to be stored. Those symlinks could look like this:

```
   Setup_2014_07_EPT_Prod -> ../../ant-database/Setup_2014_07_EPT_Prod
```

Then ../../ant-database is outside of this repository and can, for
example, be stored on blaster and sshfs'd to your local machine.

It is also a good idea to use git again in those database folders in
case something is messed up or for backup/migration purposes. However,
the database can get quite large, so it was separated from this code
repository.
