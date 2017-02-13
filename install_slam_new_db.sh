#!/bin/bash

# Copyright 2014 David Ainsworth
#
# This file is part of SLAM
#
# SLAM is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SLAM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with SLAM.  If not, see <http://www.gnu.org/licenses/>.

#set -x
#set -e

set -e
WORKING_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SLAM="$WORKING_DIR"/SLAM
if [[ $# -eq 0 ]] ; then
    echo 'Please provide directory name'
    exit 0
fi
mkdir -p $1
cd $1
if (( $# >= 1 ))
then
	shift
else
	echo "Please give install directory e.g. \"./install_slam.sh INSTALL_DIR bacteria viruses\""
	exit 1
fi

NCBI_FTP="ftp://ftp.ncbi.nih.gov/genomes/refseq"
TAX_DB_FILE_NAME="taxdump.tar.gz"
BACTERIA=false
VIRUSES=false
while (( $# > 0 ))
do
case $1 in
        "bacteria")
        BACTERIA=true
        ;;
        "viruses")
        VIRUSES=true
        ;;
        *)
        echo "Unknown option" $1
        ;;
esac
shift
done

if [ "$BACTERIA" = "true" ]
then
mkdir -p bacteria
cd bacteria
#rm -f $DB_FILE_NAME
if [ ! -e "baclibrarydownloaded" ]
then
        echo "Downloading RefSeq bacterial genomes"
        wget -O archaea_summary.txt ${NCBI_FTP}/archaea/assembly_summary.txt
	wget -O bacteria_summary.txt ${NCBI_FTP}/bacteria/assembly_summary.txt
        awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' archaea_summary.txt bacteria_summary.txt > ftpdirpaths
        awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gbff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
	cat ftpfilepaths | xargs -n 1 -P 32 curl -LO &> /dev/null
        echo "Decompressing"
	ls *.gz | xargs -n 1 -P 16 gzip -d
	rm ftpdirpaths archaea_summary.txt bacteria_summary.txt ftpfilepaths
        echo "Done"
        touch "baclibrarydownloaded"
else
        echo "Already downloaded bacterial genomes"
fi
cd ..
fi
if [ "$VIRUSES" = "true" ]
then
mkdir -p viruses
cd viruses
#rm -f $DB_FILE_NAME
if [ ! -e "virlibrarydownloaded" ]
then
        echo "Downloading RefSeq viral genomes"
        wget -O viruses_summary.txt ${NCBI_FTP}/viral/assembly_summary.txt
	awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' viruses_summary.txt > ftpdirpaths
        awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gbff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
	cat ftpfilepaths | xargs -n 1 -P 32 curl -LO &> /dev/null
        echo "Decompressing"
	ls *.gz | xargs -n 1 -P 16 gzip -d 
        rm ftpdirpaths viruses_summary.txt ftpfilepaths
        echo "Done"
        touch "virlibrarydownloaded"
else
        echo "Already downloaded viral genomes"
fi
cd ..
fi
mkdir -p taxonomy
cd  taxonomy
if [ ! -e "taxdownloaded" ]
then
        echo "Downloading RefSeq taxonomy"
        wget ${NCBI_FTP}/../../pub/taxonomy/$TAX_DB_FILE_NAME
        echo "Decompressing"
        tar xvf $TAX_DB_FILE_NAME
#        rm $TAX_DB_FILE_NAME
        echo "Done"
        touch "taxdownloaded"
else
        echo "Already downloaded taxonomy"
fi
cd ..

echo  "Creating taxonomy database"
"$SLAM" --parse-taxonomy taxonomy/names.dmp taxonomy/nodes.dmp --output-file taxDB
echo "Creating sequence database"
if [ "$VIRUSES" = "true" ]
then
if [ "$BACTERIA" = "true" ]
then
"$SLAM" --output-file database --parse-genbank bacteria/*.gbff viruses/*.gbff
else
"$SLAM" --output-file database --parse-genbank viruses/*.gbff
fi
elif [ "$BACTERIA" = "true" ]
then
"$SLAM" --output-file database --parse-genbank bacteria/*.gbff
fi

