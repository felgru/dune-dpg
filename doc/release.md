Release process
===============

To create a new relase of dune-dpg follow these steps:

1. Create a new branch for the release and set the version number in
   dune.module. Replace `$VERSION` with the version of the new release.
   ```
   git checkout -b releases/$VERSION
   sed -i "s/^Version: .*/Version: $VERSION/" dune.module
   ```
   Remember to change the version number on the master branch to
   something like $NEXTVERSION-git.
2. Check that all major changes are documented in
   [CHANGELOG.md](../CHANGELOG.md). The format of this file is based on
   [Keep a Changelog](http://keepachangelog.com/).
3. Build and test everything. Especially do a
   ```
   make build_tests
   ctest
   ```
   to see if all the unit tests still succeed.
4. Update the copyright years in COPYING.
5. Tag the release with
   ```
   git tag -s v$VERSION
   ```
   The `-s` option creates a GPG-signed annotated tag.
   You can add a small annotation text to describe the release.
   Use `git tag -n10` to see previously used annotations.
6. Push the newly created release tag and branch with
   ```
   git push origin v$VERSION releases/$VERSION
   ```
7. Finally, you can create a release tarball with
   ```
   git archive --prefix=dune-dpg-$VERSION/ -o dune-dpg-$VERSION.tar.gz v$VERSION
   ```
   Currently, we do not have a canonical place to publish release tarballs
   but such a tarball could be useful if you want to publish dune-dpg
   as supplemental material of a publication.
