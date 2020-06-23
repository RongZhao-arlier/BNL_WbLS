#!/bin/bash
grep class -B3 ${1}|grep ^\ >tmp.txt
