% Class: Template to build models which can be of type
%       (i) linear,
%       (ii) bilinear,
%       (iii) deterministic with/out additive noise, or
%       (iv) stochastic
% author: Nathalie Cauchi
% -------------------------------------------------------

classdef createModel
    % Class to build models
    
    properties
        dim;          % Model dimension
        dim_u;        % Input dimension
        dim_d;        % Noise dimension
        Type;         % Linear - 'l' / Bilinear - 'b'
        Noise;        % Deterministic -'d' / Stochastic - 's'
        A;            % A-matrix
        B;            % B-matrix
        C;            % C-matrix
        F;            % F-matrix (disturbances)
        N;            % N-matrix (bilinear)
        tr;           % Symbolic transitions function
        dt;           % Discrete time step 
        sigma;        % Variance
        dW;           % Weiner increments
    end
    
    methods
        %%%%%%%%%%% Initialisation of model %%%%%%%%%%
        function obj = InitialiseModel(obj,Type,Noise,A,B,C,F,N,tr,dt,sigma)
            
            obj.Type   = Type;
            obj.Noise  = Noise;
            obj.A      = A;
            obj.B      = B;
            obj.C      = C;
            obj.F      = F;
            obj.N      = N;
            obj.tr     = tr;
            obj.dt     = dt;
            
            if obj.Noise == 'd'
                obj.sigma  = [];
                obj.dW     = [];
            else
                obj.sigma = sigma;
                obj.dW    = [];       % Brownian increments
            end
            
            
            obj.dim    = size(A,2);
            obj.dim_u  = size(B,2);
            obj.dim_d  = size(F,2);
            
            %%% Discretise matrices
            obj.A = eye(size(A)) + obj.dt.*obj.A;
            obj.B = obj.dt.*obj.B;
            obj.F = obj.dt.*obj.F;
            if obj.Type ~= 'l'
                obj.N = obj.dt.*obj.N;
            end
            
        end
        
        %%%%%%%%%%% Creation of symbolic model %%%%%%%%%%
        function sys = createSymbModel(obj)
            sys = obj;
            
            % Defining symbols for states vector
            x='';xbar='';
            
            for i=1:obj.dim
                eval(['syms',' ','x',num2str(i),' ','real']);
                eval(['syms',' ','x',num2str(i),'bar ','real']);
                x=[x,'x',num2str(i),' '];
            end
            if x(end)==' '
                x=x(1:end-1);
            end
                     
            x_vec = eval(['[',x,']','''']);
            
            % Defining symbols for input vector
            u='';
            if obj.dim_u >0
                for i=1:obj.dim_u
                    eval(['syms',' ','u',num2str(i),' ','real'])
                    u=[u,'u',num2str(i),' '];
                end
                if u(end)==' '
                    u=u(1:end-1);
                end
            end
            u_vec = eval(['[',u,']','''']);
            
            if obj.Type == 'l'
                if obj.dim_u >0
                    sys_par = obj.A*x_vec+obj.B*u_vec;
                else
                    sys_par = obj.A*x_vec;
                end
                
            else
                if isempty(obj.N)
                    disp('Selected Bilinear model as type and N-matrix is empty!');
                    return;
                end
                sys_bi = obj.N*kron(u_vec,x_vec);
                
                sys_par = obj.A*x_vec+obj.B*u_vec + sys_bi;
                
            end
            
            % Defining symbols for disturbance vector
            d='';
            
            if obj.dim_d >0
                for i=1:obj.dim_d
                    eval(['syms',' ','d',num2str(i),' ','real'])
                    d=[d,'d',num2str(i),' '];
                end
                if d(end)==' '
                    d=d(1:end-1);
                end
                d_vec = eval(['[',d,']','''']);
                
                sys_mod       = sys_par + obj.F*d_vec;
                sys.tr=eval(['matlabFunction(sys_mod,''vars'',[',x,' ',u,' ',d,'])']);
                
                if obj.Noise == 's'
                    s='';
                    for i=1:obj.dim
                        eval(['syms',' ','s',num2str(i),' ','real']);
                        s=[s,'s',num2str(i),' '];
                    end
                    if s(end)==' '
                        s=s(1:end-1);
                    end
                    s_vec = eval(['[',s,']','''']);
                    
                    sys_mod       = sys_mod + obj.sigma*s_vec;
                    sys.tr=eval(['matlabFunction(sys_mod,''vars'',[',x,' ',u,' ',d,' ',s,'])']);
                    
                end
            else
                if obj.Noise == 's'
                    s='';
                    for i=1:obj.dim
                        eval(['syms',' ','s',num2str(i),' ','real']);
                        s=[s,'s',num2str(i),' '];
                    end
                    if s(end)==' '
                        s=s(1:end-1);
                    end
                    s_vec = eval(['[',s,']','''']);
                    
                    sys_mod       = sys_par + obj.sigma*s_vec;
                    %%%%%%%%% Finalization of the kernel %%%%%%%%%%%%
                    sys.tr=eval(['matlabFunction(sys_mod,''vars'',[',x,' ',u,' ',d,' ',s,'])']);
                else
                    sys.tr=eval(['matlabFunction(sys_par,''vars'',[',x,' ',u,'])']);
                    
                end
            end
            
            
            
        end
        %%%%%%%%% Execute model over time T %%%%%%%%%%%%
        function y = runModel(sys,x_init, U,D,  T)
            % Check if input signals are of correct dimensions
            if size(sys.A,2) ~= size(x_init,1)
                disp('Incorrect initial conditions');
                return;
            end
            
            if size(sys.B,2) ~= size(U,2)
                disp('Incorrect input signal');
                return;
            end
            if size(sys.F,2) ~= size(D)
                disp('Incorrect disturbance signal');
                return;
            end
            
            x      = zeros(size(sys.A,2),T);
            y      = zeros(size(sys.C*x,1),T);
            s      = zeros(size(sys.A,2),T);
            switch size(sys.A,2)
                case 1
                    switch size(sys.B,2)
                        case 0
                            x(:,1)=  sys.tr(x_init(1,1),D(1,1),D(1,2));
                            for i = 1:T
                                x(:,i+1) = sys.tr(x(1,i),D(i,1),D(i,2));
                                
                            end
                            y = sys.C*x;
                        case 1
                            switch size(sys.F,2)
                                case 0
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),U(1,1));
                                    else
                                        x(:,1)=  sys.tr(x_init(1,1),U(1,1),0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1));
                                        else
                                            s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),s(1,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 1
                                    if sys.Noise == 'd'
                                        x(:,1) =  sys.tr(x_init,U(1,1),D(1,1));
                                    else
                                        %sum(sys.dW(1:1));
                                        x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),sum(sys.dW(1:1)));
                                    end
                                    
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1));
                                        else
                                            s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),s(1,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 2
                                    if sys.Noise == 'd'
                                        x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2));
                                    else
                                        x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2),0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2));
                                        else
                                            s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2),s(1,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 3
                                    if sys.Noise == 'd'
                                        x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2),D(1,3));
                                    else
                                        x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2),D(1,3),0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2),D(i,3));
                                        else
                                            s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2),D(i,3),s(1,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 4
                                    if sys.Noise == 'd'
                                        x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2),D(1,3),D(1,4));
                                    else
                                        x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2),D(i,3),D(i,4));
                                        else
                                            s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2),D(i,3),D(i,4),s(1,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 5
                                    if sys.Noise == 'd'
                                        x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),D(1,5));
                                    else
                                        x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),D(1,5),0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2),D(i,3),D(i,4),D(i,5));
                                        else
                                            s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),s(1,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end                                    
                            end
                            
                            
                        case 2
                            switch size(sys.F,2)
                                case 0
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),U(1,1),U(1,2));
                                    else
                                        x(:,1)=  sys.tr(x_init(1,1),U(1,1),U(1,2),0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),U(i,2));
                                        else
                                            s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),U(i,2),s(1,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 1
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),U(1,1),U(1,2),D(1,1));
                                    else
                                        x(:,1)=  sys.tr(x_init(1,1),U(1,1),U(1,2),D(1,1),0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),U(i,2),D(i,1));
                                        else
                                            s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),U(i,2),D(i,1),s(1,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 2
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),U(1,1),U(1,2),D(1,1),D(1,2));
                                    else
                                        x(:,1)=  sys.tr(x_init(1,1),U(1,1),U(1,2),D(1,1),D(1,2),0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),U(i,2),D(i,1),D(i,2));
                                        else
                                            s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),U(i,2),D(i,1),D(i,2),s(1,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                            end
                        case 3
                            switch size(sys.F,2)
                                case 0
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),U(1,1),U(1,2),U(1,3));
                                    else
                                        x(:,1)=  sys.tr(x_init(1,1),U(1,1),U(1,2),U(1,3),0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),U(i,2),U(i,3));
                                        else
                                            s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),U(i,2),U(i,3),s(1,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 1
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),U(1,1),U(1,2),U(1,3),D(1,1));
                                    else
                                        x(:,1)=  sys.tr(x_init(1,1),U(1,1),U(1,2),U(1,3),D(1,1),0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),U(i,2),U(i,3),D(i,1));
                                        else
                                            s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),U(i,2),U(i,3),D(i,1),s(1,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 2
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),U(1,1),U(1,2),U(1,3),D(1,1),D(1,2));
                                    else
                                        x(:,1)=  sys.tr(x_init(1,1),U(1,1),U(1,2),U(1,3),D(1,1),D(1,2),0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),U(i,2),U(i,3),D(i,1),D(i,2));
                                        else
                                            s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),U(i,1),U(i,2),U(i,3),D(i,1),D(i,2),s(1,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                            end
                    end
                case 2
                    switch size(sys.B,2)
                         case 0
                             switch size(sys.F,2)
                                case 0
                                    if sys.Noise == 'd'
                                     x(:,1)=  sys.tr(x_init(1,1),x_init(2,1));
 
                                    else
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),0,0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),s(1,i),s(2,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end                          
                                case 1
                                    if sys.Noise == 'd'
                                         x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),D(1,1));
                                    else
                                         x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),D(1,1),0,0);                           
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),D(i,1));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),D(i,1),s(1,i),s(2,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                 case 2
                                    if sys.Noise == 'd'
                                         x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),D(1,1),D(1,2));
                                    else
                                         x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),D(1,1),D(1,2),0,0);                           
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),D(i,1),D(i,2));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),D(i,1),D(i,2),s(1,i),s(2,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                     
                                 case 3
                                    if sys.Noise == 'd'
                                         x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),D(1,1),D(1,2),D(1,3));
                                    else
                                         x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),D(1,1),D(1,2),D(1,3),0,0);                           
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),D(i,1),D(i,2),D(i,3));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),D(i,1),D(i,2),D(i,3),s(1,i),s(2,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                             
                             end
                        case 1
                            switch size(sys.F,2)
                                case 0
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),U(1,1));
                                    else
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),U(1,1),0,0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),s(1,i),s(2,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 1
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),U(1,1),D(1,1));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),U(1,1),D(1,1),s(1,1),s(2,1));
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),D(i,1));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),D(i,1),s(1,i),s(2,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 2
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),U(1,1),D(1,1),D(1,2));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),U(1,1),D(1,1),D(1,2),0,0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),D(i,1),D(i,2));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),D(i,1),D(i,2),s(1,i),s(2,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 3
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),U(1,1),D(1,1),D(1,2),D(1,3));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),U(1,1),D(1,1),D(1,2),D(1,3),0,0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),D(i,1),D(i,2),D(i,3));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),D(i,1),D(i,2),D(i,3),s(1,i),s(2,i));
                                        end
                                        
                                    end
                                    y = sys.C*x;   
                                case 4
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),U(1,1),D(1,1),D(1,2),D(1,3),D(1,4));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),0,0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),D(i,1),D(i,2),D(i,3),D(i,4));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),D(i,1),D(i,2),D(i,3),D(i,4),s(1,i),s(2,i));
                                        end
                                        
                                    end
                                    y = sys.C*x;
                                case 5
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),D(1,5));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),D(1,5),0,0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),D(i,1),D(i,2),D(i,3),D(i,4),D(i,5));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),s(1,i),s(2,i));
                                        end
                                        
                                    end
                                    y = sys.C*x;                                    
                            end
                            
                        case 2
                            switch size(sys.F,2)
                                case 0
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),U(1,1),U(1,2));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),U(1,1),U(1,2),0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),U(i,2));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),U(i,2),s(1,i),s(2,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 1
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),U(1,1),U(1,2),D(1,1));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),U(1,1),U(1,2),D(1,1),0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),U(i,2),D(i,1));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),U(i,2),D(i,1),s(1,i),s(2,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 2
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),U(1,1),U(1,2),D(1,1),D(1,2));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),U(1,1),U(1,2),D(1,1),D(1,2),0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),U(i,2),D(i,1),D(i,2));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),U(i,2),D(i,1),D(i,2),s(1,i),s(2,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 3
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),U(1,1),U(1,2),D(1,1),D(1,2),D(1,3));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),U(1,1),U(1,2),D(1,1),D(1,2),D(1,3),0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),U(i,2),D(i,1),D(i,2),D(i,3));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),U(i,2),D(i,1),D(i,2),D(i,3),s(1,i),s(2,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 4
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),U(1,1),U(1,2),D(1,1),D(1,2),D(1,3),D(1,4));
                                    else
                                        x(:,1)= ys.tr(x_init(1,1),x_init(2,1),U(1,1),U(1,2),D(1,1),D(1,2),D(1,3),D(1,4),0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),U(i,2),D(i,1),D(i,2),D(i,3),D(i,4));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),U(i,1),U(i,2),D(i,1),D(i,2),D(i,3),D(i,4),s(1,i),s(2,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                            end
                            
                            
                    end
                case 3
                    switch size(sys.B,2)
                        case 1
                            switch size(sys.F,2)
                                case 0
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),0, 0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i), U(i,1));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i), U(i,1),s(1,i),s(2,i),s(3,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 1
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),D(1,1));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),D(1,1),0, 0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i), U(i,1),D(i,1));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i), U(i,1),D(i,1),s(1,i),s(2,i),s(3,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 2
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),D(1,1),D(1,2));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),D(1,1),D(1,2),0, 0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i), U(i,1),D(i,1),D(i,2));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i), U(i,1),D(i,1),D(i,2),s(1,i),s(2,i),s(3,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 3
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),D(1,1),D(1,2),D(1,3));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),D(1,1),D(1,2),D(1,3),0, 0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i), U(i,1),D(i,1),D(i,2),D(i,3));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i), U(i,1),D(i,1),D(i,2),D(i,3),s(1,i),s(2,i),s(3,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 4
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),D(1,1),D(1,2),D(1,3),D(1,4));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),0, 0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i), U(i,1),D(i,1),D(i,2),D(i,3),D(i,4));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i), U(i,1),D(i,1),D(i,2),D(i,3),D(i,4),s(1,i),s(2,i),s(3,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 5
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),D(1,5));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),D(1,5),0, 0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i), U(i,1),D(i,1),D(i,2),D(i,3),D(i,4),D(i,5));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i), U(i,1),D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),s(1,i),s(2,i),s(3,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 2
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),D(1,1),D(1,2));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),D(1,1),D(1,2),0, 0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i), U(i,1),D(i,1),D(i,2));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i), U(i,1),D(i,1),D(i,2),s(1,i),s(2,i),s(3,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end                                    
                            end
                            
                            
                        case 2
                            switch size(sys.F,2)
                                case 0
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),U(1,2));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),U(1,2),0, 0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),U(i,1),U(i,2));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),U(i,1),U(i,2),s(1,i),s(2,i),s(3,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 1
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),U(1,2),D(1,1));
                                    else
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),U(1,2),D(1,1),0, 0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),U(i,1),U(i,2),D(i,1));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),U(i,1),U(i,2),D(i,1),s(1,i),s(2,i),s(3,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                                case 2
                                    if sys.Noise == 'd'
                                        x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),U(1,2),D(1,1),D(1,2));
                                    else
                                        x(:,1)= sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),U(1,1),U(1,2),D(1,1),D(1,2),0, 0, 0);
                                    end
                                    for i = 1:T
                                        if sys.Noise == 'd'
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),U(i,1),U(i,2),D(i,1),D(i,2));
                                        else
                                            s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                            x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),U(i,1),U(i,2),D(i,1),D(i,2),s(1,i),s(2,i),s(3,i));
                                        end
                                        y(:,i) = sys.C*x(:,i);
                                    end
                            end
                            
                            
                    end
                case 4
                    switch size(sys.F,2)
                        case 0 
                             if sys.Noise == 'd'
                                x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),x_init(4,1),U(1,1));
                            else
                                x(:,1)= sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),x_init(4,1),U(1,1),0, 0, 0,0);
                            end
                            
                            for i = 1:T
                                if sys.Noise == 'd'
                                    x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),x(4,i), U(i,1));
                                else
                                    s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                    x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),x(4,i),U(i,1),s(1,i),s(2,i),s(3,i),s(4,i));
                                end
                                
                                y(:,i) = sys.C*x(:,i);
                            end
                        case 1
                            if sys.Noise == 'd'
                                x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),x_init(4,1),U(1,1),D(1,1));
                            else
                                x(:,1)= sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),x_init(4,1),U(1,1),D(1,1),0, 0, 0,0);
                            end
                            
                            for i = 1:T
                                if sys.Noise == 'd'
                                    x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),x(4,i), U(i,1), D(i,1));
                                else
                                    s(:,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                    x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),x(4,i),U(i,1),D(i,1),s(1,i),s(2,i),s(3,i),s(4,i));
                                end
                                
                                y(:,i) = sys.C*x(:,i);
                            end

                        case 3
                            
                            x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),x_init(4,1),U(1,1),D(1,1),D(1,2),D(1,3));
                            
                            for i = 1:T
                                
                                x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),x(4,i), U(i,1), D(i,1),D(i,2),D(i,3));
                                
                                y(:,i) = sys.C*x(:,i);
                            end
                        case 4
                            
                            x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),x_init(4,1),U(1,1),D(1,1),D(1,2),D(1,3),D(1,4));
                            
                            for i = 1:T
                                
                                x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),x(4,i), U(i,1), D(i,1),D(i,2),D(i,3),D(i,4));
                                
                                y(:,i) = sys.C*x(:,i);
                            end
                            
                        case 5
                            
                            x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),x_init(4,1),U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),D(1,5));
                            
                            for i = 1:T
                                
                                x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),x(4,i), U(i,1), D(i,1),D(i,2),D(i,3),D(i,4),D(i,5));
                                
                                y(:,i) = sys.C*x(:,i);
                            end
                        case 6
                            
                            x(:,1)=  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),x_init(4,1),U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),D(1,5),D(1,6));
                            
                            for i = 1:T
                                
                                x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),x(4,i), U(i,1), D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6));
                            end
                            y = sys.C*x;
                    end
                    
                    
                case 7
                    switch size(sys.B,2)
                        case 1
                            if sys.Noise == 'd'
                                x(:,1) =  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),x_init(4,1),x_init(5,1), x_init(6,1), x_init(7,1),U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),D(1,5),D(1,6),D(1,7));
                            else
                                x(:,1) =  sys.tr(x_init(1,1),x_init(2,1),x_init(3,1),x_init(4,1),x_init(5,1), x_init(6,1), x_init(7,1),U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),D(1,5),D(1,6),D(1,7),0);
                            end
                            for i = 1:T
                                if sys.Noise == 'd'
                                    x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),x(4,i),x(5,i), x(6,i), x(7,i),U(i,1),D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6),D(i,7));
                                else
                                    s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                    x(:,i+1) = sys.tr(x(1,i),x(2,i),x(3,i),x(4,i),x(5,i), x(6,i), x(7,i),U(i,1),D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6),D(i,7),s(1,i));
                                end
                                y(:,i) = sys.C*x(:,i);
                            end
                        case 2
                            if sys.Noise == 'd'
                                x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2));
                            else
                                x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2),0);
                            end
                            for i = 1:T
                                if sys.Noise == 'd'
                                    x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2));
                                else
                                    s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                    x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2),s(1,i));
                                end
                                y(:,i) = sys.C*x(:,i);
                            end
                        case 3
                            if sys.Noise == 'd'
                                x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2),D(1,3));
                            else
                                x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2),D(1,3),0);
                            end
                            for i = 1:T
                                if sys.Noise == 'd'
                                    x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2),D(i,3));
                                else
                                    s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                    x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2),D(i,3),s(1,i));
                                end
                                y(:,i) = sys.C*x(:,i);
                            end
                        case 4
                            if sys.Noise == 'd'
                                x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2),D(1,3),D(1,4));
                            else
                                x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),0);
                            end
                            for i = 1:T
                                if sys.Noise == 'd'
                                    x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2),D(i,3),D(i,4));
                                else
                                    s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                    x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2),D(i,3),D(i,4),s(1,i));
                                end
                                y(:,i) = sys.C*x(:,i);
                            end
                            
                        case 5
                            if sys.Noise == 'd'
                                x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),D(1,5));
                            else
                                x(:,1) =  sys.tr(x_init,U(1,1),D(1,1),D(1,2),D(1,3),D(1,4),D(1,5),0);
                            end
                            for i = 1:T
                                if sys.Noise == 'd'
                                    x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2),D(i,3),D(i,4),D(i,5));
                                else
                                    s(1,i) = sum(sys.dW((i-1)+1:i)); % Weiner Increment
                                    x(:,i+1) = sys.tr(x(1,i),U(i,1),D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),s(1,i));
                                end
                                y(:,i) = sys.C*x(:,i);
                            end
                            
                    end
            end
            
            
            
        end
        
        function S = stepResponse(sys,U,D,T)
            x      = zeros(size(sys.A,2),T);
            nu     = size(U,2);
            x_init = x(:,1);
            
            
            for i = 1:nu
                x(:,1) = zeros(size(sys.A,2),1);
                U(:,i) = ones(T,1);
                U(1:5,i) = zeros(5,1);
                y = runModel(sys, x_init,U,D,T);
                plot(1:T,y);
                if size(y,1) == 1
                    S(i) = stepinfo(y,1:T,y(:,end));
                else
                    S(i,1) = stepinfo(y(1,:),1:T,y(1,end));
                    S(i,2) = stepinfo(y(2,:),1:T,y(2,end));
                end
                
            end
            
        end
        
        
        
    end
    
end
