#!/bin/sh
 
# git add all files that have been modified by doctoc in the pre-commit hook
# and commit them 

PATH=$PATH:/usr/local/bin:/usr/local/sbin

# Go to git root directory
rootdir=`git rev-parse --show-toplevel`
cd $rootdir

echo
if [ -a .commit ]
then
  modified_by_doctoc=`cat .commit`
  echo "Updated ToC for:"
  echo $modified_by_doctoc
  rm .commit
  git add $modified_by_doctoc
  git commit --amend -C HEAD --no-verify --allow-empty
fi
exit
