function sFunModell(block)

	setup(block);

end


% *************************************************************************
% Initialisierung
% *************************************************************************
function setup(block)

	% Anzahl der Ein-/Ausg�nge
	block.NumInputPorts  = 1;
	block.NumOutputPorts = 2;
% Test git

	% Eigenschaften des Eingangs
	block.InputPort(1).Dimensions	= 1;
	block.InputPort(1).DatatypeID	= 0;  % double
	block.InputPort(1).Complexity	= 'Real';
	block.InputPort(1).DirectFeedthrough = false;
	block.InputPort(1).SamplingMode = 'Sample';

	% Eigenschaften des 1. Ausgangs
	block.OutputPort(1).Dimensions	= 1;
	block.OutputPort(1).DatatypeID	= 0; % double
	block.OutputPort(1).Complexity	= 'Real';
	block.OutputPort(1).SamplingMode = 'Sample';

	% Eigenschaften des 2. Ausgangs
		block.OutputPort(2).Dimensions	= 1;
	block.OutputPort(2).DatatypeID	= 0; % double
	block.OutputPort(2).Complexity	= 'Real';
	block.OutputPort(2).SamplingMode = 'Sample';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
	% Anzahl der Zust�nde
	block.NumContStates =4 ;

	% Anzahl der Parameter
	block.NumDialogPrms =7 ;


	% Abtastzeit definieren -> zeitkontinuierlich
	block.SampleTimes = [0 0];


	
	% weitere Methoden registrieren
	block.RegBlockMethod('InitializeConditions', @InitializeConditions);
	block.RegBlockMethod('Outputs', @Outputs);
	block.RegBlockMethod('Derivatives', @Derivatives);
	block.RegBlockMethod('Terminate', @Terminate);

end


% *************************************************************************
% Anfangsbedingungen setzen
% *************************************************************************
function InitializeConditions(block)

	block.ContStates.Data = [pi/1.01 0 pi/1.01 0];
	
end


% *************************************************************************
% Ausg�nge berechnen
% *************************************************************************
function Outputs(block)

	% Zust�nde auslesen
	x		= block.ContStates.Data;
    phi1=x(1);
    dphi1=x(2);
    phi2=x(3);
    dphi2=x(4);

	block.OutputPort(1).Data = phi1;
	block.OutputPort(2).Data = phi2;

end


% *************************************************************************
% Ableitungen berechnen
% *************************************************************************
function Derivatives(block)

	% Parameter auslesen
	g = block.DialogPrm(1).Data;
	m1 = block.DialogPrm(2).Data;
	m2 = block.DialogPrm(3).Data;
	l1 = block.DialogPrm(4).Data;
	l2 = block.DialogPrm(5).Data;
	Rp1 = block.DialogPrm(6).Data;
    Rp2 = block.DialogPrm(7).Data;

	% Zust�nde auslesen
	x		= block.ContStates.Data;
	phi1	= x(1);
	dphi1   = x(2);
	phi2    = x(3);
	dphi2	= x(4);

	% Eingang auslesen
	M = block.InputPort(1).Data(1);

	
	
	% Ableitungen berechnen
	ddphi1 = -(6*(4*Rp1*dphi1*l2 - 4*M*l2 + 4*Rp2*dphi1*l2 - 4*Rp2*dphi2*l2 + 6*Rp2*dphi1*l1*cos(phi1 - phi2) - 6*Rp2*dphi2*l1*cos(phi1 - phi2) + (3*g*l1*l2*m2*sin(phi1 - 2*phi2))/2 + (3*dphi1^2*l1^2*l2*m2*sin(2*phi1 - 2*phi2))/2 + 2*dphi2^2*l1*l2^2*m2*sin(phi1 - phi2) + 2*g*l1*l2*m1*sin(phi1) + (5*g*l1*l2*m2*sin(phi1))/2))/(l1^2*l2*(8*m1 + 15*m2 - 9*m2*cos(2*phi1 - 2*phi2))) ;
    
    
    
	ddphi2 =  (6*(4*Rp2*dphi1*l1*m1 - 6*M*l2*m2*cos(phi1 - phi2) + 12*Rp2*dphi1*l1*m2 - 4*Rp2*dphi2*l1*m1 - 12*Rp2*dphi2*l1*m2 + 6*dphi1^2*l1^2*l2*m2^2*sin(phi1 - phi2) - 3*g*l1*l2*m2^2*sin(phi2) + (3*dphi2^2*l1*l2^2*m2^2*sin(2*phi1 - 2*phi2))/2 + 3*g*l1*l2*m2^2*sin(2*phi1 - phi2) + 6*Rp1*dphi1*l2*m2*cos(phi1 - phi2) + 6*Rp2*dphi1*l2*m2*cos(phi1 - phi2) - 6*Rp2*dphi2*l2*m2*cos(phi1 - phi2) + 2*dphi1^2*l1^2*l2*m1*m2*sin(phi1 - phi2) - (g*l1*l2*m1*m2*sin(phi2))/2 + (3*g*l1*l2*m1*m2*sin(2*phi1 - phi2))/2))/(l1*l2^2*m2*(8*m1 + 15*m2 - 9*m2*cos(2*phi1 - 2*phi2)));
	
	

	% Ableitungen zuweisen
	block.Derivatives.Data = [dphi1; ddphi1; dphi2; ddphi2];
	
end


% *************************************************************************
% Aufr�umen (wenn n�tig)
% *************************************************************************
function Terminate(block)
end
