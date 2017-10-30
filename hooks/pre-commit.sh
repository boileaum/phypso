#!/bin/sh

PATH=$PATH:/usr/local/bin:/usr/local/sbin

# Go to git root directory
rootdir=`git rev-parse --show-toplevel`
cd $rootdir

# List all .md files that do not contain a DOCTOC SKIP comment and apply doctoc
modified_by_doctoc=`git diff --cached --name-only |ack -xL 'DOCTOC SKIP' | grep '\.md$'`

rm -f .commit
# Test if list is not empty (remove spaces)
if [ -n "${modified_by_doctoc// }" ]
then
  echo "pre-commit hook is running doctoc"
  doctoc --notitle --gitlab $modified_by_doctoc
  # Store list in a tmp file for post-commit hook
  echo $modified_by_doctoc > .commit
fi

