#snipSNPer.py README#

Open terminal and execute $ python snipSNPer.py to launch program

Requirements:
1. Python 2.6 or later (but not Python 3.x)
	Libraries:
	1. rpy2
	2. Biopython
2. R
	Packages:
	1. ggplot2
3. primer3
4. Internet access

##Installation:

1. To install the required Python packages, use the pip installation tool by calling "$ sudo pip install [package name]" from the command line
2. Install R with the ggplot2 package using install.pachages("ggplot2") from inside an R session
3. Install primer 3 by downloading the latest version from http://primer3.sourceforge.net/releases.php
	1. Follow installation instructions for **LINUX/UNIX, then MAC** on http://primer3.sourceforge.net/primer3_manual.htm#installLinux
	**IT IS VERY IMPORTANT TO FOLLOW THE LINUX/UNIX INSTRUCTIONS FIRST, THEN THE MAC INSTRUCTIONS, OTHERWISE THE FILES NEEDED IN THE MAC INSTRUCTIONS WILL NOT EXIST**
	2. Copy the primer3_config folder in the downloaded primer3 folder into the /opt folder on the computer
	3. To test, change directories into the downloaded primer3 folder and execute the command $ primer3_core -output=out.txt example
	4. The out.txt file should contain no errors

EVERY PROMPT IS CASE AND FORMAT SENSITIVE. IF SOMETHING ISN'T WORKING, MAKE SURE THE CASE AND FORMAT IS CORRECT.

##Instructions:

1. When prompted, enter the the genomic point or region around which you are looking for snipSNP sites. Points and regions must be entered in the form chromosome:position (ex. I:1000000 for basepair 1000000 on chromosome I or I:1000000..1100000 for the region from 1000000 to 1100000 on chromosome 1).

2. When prompted, enter the way in which you would like to enter the enzymes that you want to use. "Standards" draws the TimeSaver enzymes from NEB. "List" will prompt you to enter the enzymes separated by commas. "File" will prompt you to select a .csv file with the desired enzymes.

3. When prompted, enter the two strains that you would like to use. Order is not important.

4. Saving the chromosome map:
	1. When prompted, enter whether you would like to save the chromosome map. If you don not want to save it, the program will assign it a temporary name and remove it at the end of the program.
	2. When prompted, enter a file name.
	3. File will open automatically with Preview

5. Enter the numbers of the sites that you would like to use. If there are 5 sites on either side of the point/region, sites 5 and 6 will be the two immediately outside the point/region.

6. When prompted, enter a name for the output .csv file.

7. When you see: "Final file output to" and your filename, the program is done running.