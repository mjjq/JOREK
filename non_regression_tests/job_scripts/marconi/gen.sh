#!/bin/bash
for i in kink_*job tear_*job; do NAME=${i/.job/}; sed  "s/_CASE_/$NAME/g" case.skel > $i; done
chmod +x *.job
