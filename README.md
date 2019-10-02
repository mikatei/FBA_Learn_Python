# FBA Learn Python


This module allows you to run and make examples of FBA.


To install required Python packages:

```
$bash install_dependencies.sh
```

To run the program:

```
$python main.py
```


The 4 required inputs to the program are:  
*(1): a list of reactions as exemplified by the files in the 'Examples' directory. 

*(2): An index for the variable you want to optimize (objective function) 

*(3): The type of optimization (min or max) 

*(4): The max flux value through a reaction (e.g. 100 or 1000) 

The output will be:
*(1) The stoichiometric matrix made out of the list of reactions.
*(2) The flux vector with all the flux values. 



Make an FBA reaction model and put the file in the 'Examples' directory.
* The model should look like all the other reaction models in the Examples directory (take a look), with the reaction name (or number) on the left (make sure not to repeat reaction names), then a colon ":", and then the reaction, with each compound preceded by its stoichiometric coefficient (e.g "1", "-1", or whatever integer it may be, "-37?").
A requirement for naming the file is that it doesn't repeat another file name.
Make sure to remember the exact name of the file. For example, you can name it "myexample1.txt"

Then to run the program run 

```
"python main.py"
```

The user prompt will ask you for the name of your file.


Then the user prompt will ask you for the index of the reaction whose flux you want to maximize.
( e.g. If the reaction is the third one you listed in the file, write "3")

Continue following the user prompt and you'll get the solved flux vector.

In order to run the program quickly without the user input, go into main.py, and in the function 'main():', comment out the user input function,
and uncomment the quick process function, replace "myexample1.txt" with your file name, replace all the variables with yours, and run "main.py" again.










