# FBA Learn Python


This module allows you to run and make examples of FBA.
It uses 'optlang' to solve the linear programming optimization problems.


Once you've downloaded the files from github, you can run ". install_dependencies.sh" to install the necessary modules to run the program.



Definitions:
The format of our reactions is shown in the 'Examples' directory. Every text file there is a list of reactions in a specific format.
When you add your Model, please make it in the same format- otherwise the program will not run properly.

How to Run The Program:
Make an FBA reaction model and put the file in the 'Examples' directory.
* The model should look like all the other reaction models in the Examples directory (take a look), with the reaction name (or number) on the left (make sure not to repeat reaction names), then a colon ":", and then the reaction, with each compound preceded by its stoichiometric coefficient (e.g "1", "-1", or whatever integer it may be, "-37?").
A requirement for naming the file is that it doesn't repeat another file name.
Make sure to remember the exact name of the file. For example, you can name it "myexample1.txt"

Then to run the program run "python main.py"

The user prompt will ask you for the name of your file.

Then the user prompt will ask you for the index of the reaction whose flux you want to maximize.
( e.g. If the reaction is the third one you listed in the file, write "3")

Continue following the user prompt and you will get the solved flux vector.

In order to run the program quickly without the user input, go into main.py, and in the function 'main():', comment out the user input function,
and uncomment the quick process function, replace "myexample1.txt" with your file name, replace all the variables with yours, and run "python main.py" again.









In order to run the program from inside the code:

You can make an FBA reaction model and put the file in the 'Examples' directory.
* The model should look like all the other reaction models in the Examples directory (take a look), with the reaction name (or number) on the left (make sure not to repeat reaction names), then a colon ":", and then the reaction, with each compound preceded by its stoichiometric coefficient (e.g "1" or "-1").

Then go into the file Aux/aux_2.py, and at the last function of the file (the bottom), make a new variable with the location of your file and add it
to the list of filepaths, add it to the end of the list.

Then you can choose your objective variable in main- if you want the ith reaction to be the one with maximized flux, under objective_index, punch in "i-1"



