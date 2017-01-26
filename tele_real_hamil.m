%
% 18.01.2017
%-------------
% Teleportation of a single qubit state 
% a |0> + b|1> over an array of N qubits 
% by just applying a sequence of pi pulses
%
% Hamiltonian of the system
% H = sum_k Omega_k * sigma_x^(k) + sum_k Delta_k * n_k + sum_k,m V_{km} *
% n_k * n_m/|rk-rm|^gamma

% Omega_k : Rabi frequency of particle k
% Delta_k : Detuning of the rydberg state of particle k
% V_{km}  : interaction potential
% V_km = V_0 * 1/|k-m|^gamma
%Approximation
% Delta_k = 0 -> no Detuning
% H = sum_k Omega_k * sigma_x^(k) + sum_k,m V_{km} *
% n_k * n_m/|rk-rm|^gamma

clear all;
format short;



%N   = 7; % total number of atoms


% Parameters  
V_0   = 20; % nearest neighbot interaction strength
gamma = 6;  %  exponent of interaction force


% Single particle operators and states

v0  = [0;1];
v1  = [1;0];


% Pauli operators
sigma_x = [0,1;1,0];
sigma_y = [0,i;-i,0];
sigma_z = [-1,0;0,1];

% Initial state
b = 1;
a = 1;
Psi = [b;a];
Psi = Psi/norm(Psi);

one     = eye(2);   

% 3 particle operators
sig_2 = kron (sigma_y,one); 
sig_2 = sparse(kron(one,sig_2)); % 3 particle sigma operator 
                                 % 1 x sigma x 1
                                 
% number operator                                 
num     = [1,0; 0 ,0];
                              
num_1 = kron(one,one);
num_1 = sparse(kron(num,num_1));


num_3 = kron(one,num);
num_3 = sparse(kron(one,num_3));

one_N = kron(one,one);
one_N = sparse(kron(one_N,one));


sigma_x_N = kron(one,sigma_x);

Omega_vec=[0.001:0.01:1 1.1:0.1:20]; % Different Rabi frequencies
fid(1:6,1:length(Omega_vec)) = 0; % Fidelity of the teleportation process

integer = 1;

for N = 4:8% nr of partivles
    
        a = 1.0; % spacing between the spins (same disctance between the two chains)
        x = a*(1:N); % spin chain

        init = v0;
        for cnt = 1:N-2
             init  = sparse((kron(init,v0))); % initial state vector
             sigma_x_N = sparse(kron(one,sigma_x_N));
        end
        init = sparse(kron(Psi,init));

        % Sigma operators
        G_c2 = sparse(kron(sigma_y,one) ); 
        % 
        G_c1 = sparse(kron(one,sigma_y) );
        G_3 = sparse( sig_2 ); % acting on atom j-1,j,j+1


        %  Construction of number operators 
        num_op_m = sparse(N*2^N,2^N); % all N-particle number operators n_k are stored in this matrix
        help_vec_n(1:2,1:2*N) = 0;    % helping vector to produce the N particle sigma_x operator

        for cnt = 1:N 

              for cnt2 = 1:N

                    if cnt2 == cnt
                             help_vec_n(:,2*cnt2-1:2*cnt2) = num;
                             continue
                    end
                    help_vec_n(:,2*cnt2-1:2*cnt2) = one;

              end

              n_k = sparse(help_vec_n(:,2*N-1:2*N));
              for cnt3 = N-1:-1:1

                         n_k =sparse(kron(help_vec_n(:,2*cnt3-1:2*cnt3),n_k));
              end

              num_op_m((cnt-1)*2^N+1:cnt*2^N,:) = n_k;


        end

        V_int_o = sparse(N*2^N,2^N); %interaction Hamiltonian acting on the odd atoms

        for cnt = 1:N
                     V_int0 = sparse(2^N,2^N);

                     for cnt2 = 1:N

                              if cnt2 == cnt
                                   continue
                              end

                              V_int0  = sparse( V_int0 + V_0/(abs(x(cnt)-x(cnt2))^gamma) * ...
                                     num_op_m((cnt-1)*2^N+1:cnt*2^N,:)*num_op_m((cnt2-1)*2^N+1:cnt2*2^N,:));


                     end
                              V_int_o((cnt-1)*2^N+1:cnt*2^N,:) = 1/2 *V_int0;

        end

        for counter = 1:length(Omega_vec)

            Omega = Omega_vec(counter);
            t = pi/(2*Omega);


            if N > 3

                for cnt =1:N-1 % determination of Toffoli gates





                    Gate = one;
                    G_2 = one;



                    if cnt == N-1  % last block Toffoli

                        Gate = G_3;
                        G_2 = G_c1;
                        for cnt2 = 1:N-3
                            Gate = sparse(kron(one,Gate));

                            G_2 = sparse(kron(one,G_2));

                        end

                        G_2 = sparse(kron(one,G_2));


                        Psif = sparse(expm(-i*t*(Omega*G_2 +V_int_o((N-1)*2^N+1:N*2^N,:)) )*Psif);
                        Psif = sparse(expm(-i*t*( Omega*GateN +V_int_o((N-2)*2^N+1:(N-1)*2^N,:)))*Psif);

            %             
            % 
            %                 
                   elseif cnt == 1  % 1st block of Toffoli (correct)

                        for cnt2 = 1:N-4
                            Gate = sparse( kron(Gate,one));
                            G_2 = sparse(kron(G_2,one));
                        end
                        Gate = sparse(kron(G_3,Gate)); % 1st Toffoli

                        G_2 = sparse(kron(G_2,one)); %1st CNOT
                        G_2 = sparse(kron(G_c2,G_2)); 


                        Psif = sparse(expm(-i*t*(Omega*Gate + V_int_o(2^N+1:2*2^N,:)))*init); %Toffoli
                        Psif = sparse(expm(-i*t*(Omega*G_2 + V_int_o(1:2^N,:)))*Psif); %CNOT
                        GateN = Gate;



                    elseif cnt == N-2 % 2nd last block of Toffoli

                        Gate = G_3;
                        for cnt2 = 1:N-3
                            Gate = sparse(kron(one,Gate));
                        end

                         Psif = sparse(expm(-i*t*(Omega*Gate +  V_int_o((N-2)*2^N+1:(N-1)*2^N,:)) )*Psif);
                         Psif = sparse(expm(-i*t*(Omega*GateN +  V_int_o((N-3)*2^N+1:(N-2)*2^N,:)) )*Psif);
                         GateN = Gate;

                    else % cnt = 1:N-1


                     for cnt2 = N-3:-1:1 %needed for N>4

                         if cnt2 == cnt  %FILL IN THE MISSING PART


                             Gate = sparse(kron(G_3,Gate));

                             continue
                         end
                         Gate = sparse(kron(one,Gate));


                     end
                         Psif = sparse(expm(-i*t*(Omega*Gate + V_int_o(cnt*2^N+1:(cnt+1)*2^N,:)) )*Psif);
                         Psif = sparse(expm(-i*t*(Omega*GateN +  V_int_o((cnt-1)*2^N+1:(cnt)*2^N,:)) )*Psif);

                     GateN = Gate;




                    end

                end
            end


        M00 = one;
        for cnt = 1:N-1
            M00 = sparse(kron(v0,M00));
        end

        Psif = (M00' *Psif);


          fid(integer,counter) = abs(Psi'*(sigma_z)^(N-1)*((sigma_x)^(N-1))*Psif);

        end
        integer = integer +1;
end
figure()
 plot(Omega_vec/V_0,fid(1,:),'r',Omega_vec/V_0,fid(2,:),'g',Omega_vec/V_0,fid(3,:),'c',...
    Omega_vec/V_0,fid(4,:),'b', Omega_vec/V_0,fid(5,:),'k')
grid on;
grid minor;
legend('L=4 ','L=5','L=6','L=7','L=8')
xlabel('\Omega/V_0')
ylabel('Fidelity F')


