clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
R = 110;
v= -20;

% define the target's initial position and velocity. Note : Velocity
% remains contant
 


%% FMCW Waveform Generation



Rmax = 200; %m
res = 1; % m
maxV = 100; % m/s
c = 3*10^8; % Speed of light



%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

Bandwith = c/(2*res);
Tchirp = 5.5*(2*Rmax)/c;
Slope = Bandwith/Tchirp;



%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq (hz)


%a = Bandwith/Tchirp;
  
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t) % array Start from         
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    
    if i > 1
       R = R + v*(t(i) - t(i-1) ); 
    end
        
    td = 2*R/c;
    
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i)   = cos(2*pi*(fc*t(i)+(Slope*t(i)^2)/2));
    Rx (i)  = cos(2*pi*(fc*(t(i)-td)+(Slope*(t(i)-td))^2)/2);
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = cos(2*pi*(Slope*td*t(i)+(2*fc*v/c*t(i)) ));

    
end

%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix, Nr, Nd);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
signal_fft = fft(Mix(:,1));

 % *%TODO* :
% Take the absolute value of FFT output
signal_fft = abs(signal_fft/Nr);
 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
P1  = signal_fft(1:Nr/2+1) ; 

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

 % *%TODO* :
 % plot FFT output 
plot(P1) 
 
axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 10;
Td = 8;
% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 4;
Gd = 4;
% *%TODO* :
% offset the threshold by SNR value in dB
offset = 17;
% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(Nr/2,Nd);


% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR
   
   
 CFAR = RDM;  
   
for i = Tr+Gr+1 :Nr/2-(Tr+Gr)
    for j = Td+Gd+1:Nd-(Td+Gd)
        CUT = RDM(i,j);
        k = 0;
        trashold = 0;
        for x = i-Tr-Gr: i+Tr+Gr
            for y = j-Td-Gd: j+Td+Gd
               
               if (x<(i-Gr)) || (x>(Gr+i)) || (y<(i-Gd)) || (y>(Gd+i)) 
                   k = k +1;
                   trashold =  trashold + db2pow(RDM(x,y));
               end
            end
        end
        
        trashold = pow2db(trashold/k)+ offset;
        noise_level(i,j) = trashold;
        if (CUT > trashold )
            CFAR(i,j) = 1; 
        else
            CFAR(i,j) = 0;
        end
    end
end

   





% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 


for i = 1:Nr/2
    for j =1:Nd
        if (i < (Tr+Gr+1) || i > (Nr/2-(Tr+Gr))  || j < (Td+Gd+1) || j > (Nd-(Td+Gd)) )
            CFAR(i,j) = 0;
        end
    end
end





% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,CFAR);
colorbar;


 
 