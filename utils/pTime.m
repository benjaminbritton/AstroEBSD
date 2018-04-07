function text_out=pTime(text1,t1)
%simple timer code to report the time since the start
%uses the system clock function
%
%INPUT - 
%text1= string with text to display
%t1 = date vector as output from clock

%get current time
t2=clock;

%calculate the difference in time, in s, between clocks
times=floor(etime(t2,t1));

%convert to m, h and s
timem=floor(times/60);
timeh=floor(timem/60);
timed=floor(timeh/24);

%subtract off the h, m, and s
timeh2=timeh-timed*24;
timem2=timem-timed*24*60-timeh2*60;
times2=times-timed*24*60*60-timeh2*60*60-timem2*60;

%output the text
text_out=['Time since start =  [' sprintf('%3.0f',timed) ' d ' sprintf('%2.0f',timeh2) ' h ' sprintf('%2.0f',timem2) ' m ' sprintf('%2.0f',times2)  's] - ' text1];

disp(text_out);



