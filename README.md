# Barcode To Jag Files Conversion Algorithm

1. Run homology and get life and death values.
<br>
 dimension     Birth      Death
      <br>
      0 0.0000000 0.30000000
      <br>
      0 0.0000000 0.12601158
      <br>
      0 0.0000000 0.12454177
      <br>
      .....
      <br>
      1 0.2668731 0.28484301

2.    Get epsilon interval/step value. Start with 0 and then add that to take steps.
3. Check to see if anything is dead at the end of the interval's step  (Use the format to your benefit and round each death and birth value to the epsilon's decimal digits and sort them and output them accordingly)=> inclusive
	* For example, if a death is 0.032 and the epsilon is 0.01 then the death occurs during the 0.04 interval.
4.  When you note a death also note its birth and put it into the interval system developed in  step 3
5. Also note if they are betti 1 or betti 0
6. Generate the barcodes(make sure the path and naming conventions are followed)
7. Use the bar codes to generate the jag files:
	* Divided into two phases -> one for betti 1 and the other for betti 0
	* Take one barcode
	* Read it in
	* Look at the specific phases points (betti 1 will always be at the bottom and betti 0 will be at top)
	* For betti 0, all of them are alive at the start so add up all the betti 0 points -> first value in the output file's line for a specific patient and segment
	* Next go interval by interval and add up the points(subtract from a counter that holds all total number of alive betti 0 points) and write these values for each interval(go through each interval incrementally 
	* For betti 1, incrementally go through intervals and note down any born loops and add them to a counter
	* right after counting up alive loops, if any deaths are seen then then decrement the counter
	* write this counter to the line in the output file(make sure the path and naming conventions are followed)
