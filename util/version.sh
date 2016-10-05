#!/bin/bash
export LANG=C
export LC_ALL=C

USER=$PWD
RCS_VERSION=""

#RCS_VERSION="$(svnversion -n . 2>/dev/null)"
RCS_VERSION="$(svn info | grep "Revision:" | sed -e 's/^Revision: *//')"
if echo "$RCS_VERSION" 2>/dev/null | grep -i "Unversioned" 2>/dev/null >/dev/null
then RCS_VERSION=""
fi

if echo "$RCS_VERSION" 2>/dev/null | grep -i "exported" 2>/dev/null >/dev/null
then RCS_VERSION=""
fi

# True, if <STRING> is not empty
if [ -n "$RCS_VERSION" ]
then RCS_VERSION="SVN revision        ${RCS_VERSION}"
fi

# True, if <STRING> is empty
if [ -z "$RCS_VERSION" ]
then
    RCS_VERSION="$(git describe --always --dirty --abbrev 2> /dev/null)"
    if [ -n "$RCS_VERSION" ]
    then RCS_VERSION=" GIT revision        ${RCS_VERSION}"
    fi
fi

if [ -z "$RCS_VERSION" ]
then RCS_VERSION="Not under SVN or GIT version control"
fi

echo "#define RCS_VERSION           '${RCS_VERSION}'"
echo "#define COMPILE_DATE          '$(date '+%F %T')'"
