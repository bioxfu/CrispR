# CrispR Workflow Tutorial

### 1. Make project directory
```
# the project directory contains specific GROUP name and current DATE
GROUP=XXX
DATE=`date +"%Y%m%d"`
mkdir ~/Project/${GROUP}_${DATE}
```

### 2. Clone the repository and init the project
```
cd ~/Project/${GROUP}_${DATE}
git clone https://github.com/bioxfu/CrispR
cd CrispR
. init.sh
```

### 3. Copy/Download the raw data
```
DATAPATH=/the/path/of/the/raw/data/on/HPC

# if you are working on the HPC, copy the raw data
./script/copy_rawdata.sh $DATAPATH

# if you are working on the local machine, download the raw data
./script/copy_rawdata.sh $DATAPATH --download

# manually reorganize the raw data, for example:
# ./fastq/projectID/sampleID/*.gz
# each sampleID folder contains two (R1 and R2) gzipped fastq file
```

### 4. Rename the raw data
```
# dry run to check if mv command is correct
./script/rename_rawdata.sh --dry-run

# then do it 
./script/rename_rawdata.sh
```

### 5. Upload SampleSheet files to *tables* folder
```
# see example/sample_sheet.xlsx for details
# the 3rd sheet in the Excel file should be the BARCODE information
```

### 6. Upload gRNA.fa to *guide* folder 
```
# see example/gRNA.fa for details
```

### 7. Create *workflow.sh* and set the configurations in the files
```
## Arabidopsis thaliana
cp example/workflow.sh workflow_ath.sh

## Solanum lycopersicum
cp example/workflow.sh workflow_sly.sh

## edit workflow.sh
```

### 8. Run *workflow.sh*
```
# Don't run it on head node
nohup ./workflow_ath.sh &
nohup ./workflow_sly.sh &
```

### Note
```
# I update the CrispVariants to call the SNVs which are 23nt upstream and 6nt downstream from the cut site. 
# To build the CrispVariants
$MYHOME/R/$RVERSION/bin/R CMD build --no-build-vignettes CrispRVariants/
# To INSTALL the CrispVariants
$MYHOME/R/$RVERSION/bin/R CMD INSTALL -l $MYHOME/R/$RVERSION/lib64/R/library CrispRVariants_1.8.0.tar.gz 
```