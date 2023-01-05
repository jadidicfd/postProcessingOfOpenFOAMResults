delimiterIn = ' ';   % data MUST be separeted form each other by " space"
headerlinesIn = 15;  % variable start to record after line 99. the lines 
                     %before 99 are location of probes
myStructure = importdata('p1',delimiterIn,headerlinesIn);
data_p1 = myStructure.data;
myStructure = importdata('T1',delimiterIn,headerlinesIn);
data_T1 = myStructure.data;
myStructure = importdata('U1',delimiterIn,headerlinesIn);
data_U1 = myStructure.data;

myStructure = importdata('p2',delimiterIn,headerlinesIn);
data_p2 = myStructure.data;
myStructure = importdata('T2',delimiterIn,headerlinesIn);
data_T2 = myStructure.data;
myStructure = importdata('U2',delimiterIn,headerlinesIn);
data_U2 = myStructure.data;

data_p = [data_p1;data_p2];
data_T = [data_T1;data_T2];
data_U = [data_U1;data_U2];