#!/usr/bin/env bash

for i in {6..20}
do
	#./Main.py part$i > part_$i.out &
	./TRIPS_extraction.py part$i > TRIPS_part_$i &
done

