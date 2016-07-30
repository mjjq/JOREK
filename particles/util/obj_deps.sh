#!/usr/bin/env bash
make --dry-run --always-make --debug=b $1 |\
	grep 'Must remake target' |\
	grep -o "'.*'" |\
	sed -e "s/^'.mod/.obj/" -e "s/mod'/o'/" |\
	grep '[.]' |\
	tr -d "'" 
