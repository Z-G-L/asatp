#!/bin/sh


perl ../runASATP.pl --gtf anno.gtf  --trExpFile isoform_exp.txt  --output output_svg --graph yes --graphFormat svg

perl ../runASATP.pl --gtf anno.gtf  --trExpFile isoform_exp.txt  --output output_png --graph yes --graphFormat png
