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

set -e
WORKING_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SLAM="$WORKING_DIR"/SLAM
mkdir -p $1
cd $1
if (( $# >= 1 ))
then
	shift
else
	echo "Please give install directory e.g. \"./install_slam.sh INSTALL_DIR bacteria viruses\""
	exit 1
fi
NCBI_FTP="ftp://ftp.ncbi.nih.gov"
DB_FILE_NAME="all.gbk.tar.gz"
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
mkdir -p taxonomy
cd  taxonomy
if [ ! -e "taxdownloaded" ]
then
        echo "Downloading RefSeq taxonomy"
        wget ${NCBI_FTP}/pub/taxonomy/$TAX_DB_FILE_NAME
        echo "Decompressing"
        tar xvf $TAX_DB_FILE_NAME
#        rm $TAX_DB_FILE_NAME
        echo "Done"
        touch "taxdownloaded"
else
        echo "Already downloaded taxonomy"
fi
cd ..

if [ "$BACTERIA" = "true" ]
then
mkdir -p bacteria
cd bacteria
#rm -f $DB_FILE_NAME
if [ ! -e "baclibrarydownloaded" ]
then
	echo "Downloading RefSeq bacterial genomes"
	wget ${NCBI_FTP}/genomes/Bacteria/$DB_FILE_NAME
	echo "Decompressing"
	tar xf $DB_FILE_NAME
#	rm $DB_FILE_NAME
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
        wget ${NCBI_FTP}/genomes/Viruses/$DB_FILE_NAME
        echo "Decompressing"
        tar xf $DB_FILE_NAME
#       rm $DB_FILE_NAME
        echo "Done"
        touch "virlibrarydownloaded"
else
        echo "Already downloaded viral genomes"
fi
cd ..
fi

echo  "Creating taxonomy database"
"$SLAM" --parse-taxonomy taxonomy/names.dmp taxonomy/nodes.dmp --output-file taxDB
echo "Creating sequence database"
if [ "$VIRUSES" = "true" ]
then
if [ "$BACTERIA" = "true" ]
then
"$SLAM" --output-file database --parse-genbank bacteria/*/*.gbk viruses/*/*.gbk
else
"$SLAM" --output-file database --parse-genbank viruses/*/*.gbk
fi
elif [ "$BACTERIA" = "true" ]
then
"$SLAM" --output-file database --parse-genbank bacteria/*/*.gbk
fi
#rm -rf taxonomy bacteria
