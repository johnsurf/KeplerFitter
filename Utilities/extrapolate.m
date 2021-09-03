        function [Rpos Rdot, Jacobian, Hessian] = extrapolate(this, t)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            DU    = 1.0;
            VU    = 1.0;
            TU    = 1.0;
            % if nargin >= 3
            %     TU = this.TU;
            %     DU = this.DU;
            %     VU = this.VU;
            %     AU = this.AU;
            % end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            % e          = TrackParams(1);
            % p          = TrackParams(2);
            % I          = TrackParams(3);
            % omega      = TrackParams(4);
            % Omega      = TrackParams(5);
            % Mp         = TrackParams(6);
            % a          = TrackParams(7);
            % n          = TrackParams(8);
            % Pvec       = TrackParams(9:11);
            % Qvec       = TrackParams(12:14);
            % %Wvec      = TrackParams(15:17);
            % %minT       = TrackParams(18);
            % %maxT       = TrackParams(19);

            e          = this.e;
            p          = this.p;
            I          = this.Inclination;
            omega      = this.omega;
            Omega      = this.Omega;
            Mp         = this.Mp;
            a          = this.a;
            n          = this.meanMotion;
            Pvec       = this.Pvec;
            Qvec       = this.Qvec;
            %Wvec      = this.Wvec;
            %minT       = this.minT;
            %maxT       = this.maxT;
            
            onePe = 1.0 + e;                    % 39)
            oneMe = 1.0 - e;                    % 40)
            fac   = onePe*oneMe;                % 41)
            if (fac > 0.0)
                rootfac = sqrt(fac);            % 42)
            else
                rootfac = 0;                    %
            end
            
            % invert Kepler's Equation
            % Brute Force iteration (good way to seed Newton's method which follows)
            % M = MeanAnomalyEpoch;
            %  meanAnomaly M to eccentric Anomaly E
            M = n*(t/TU) + Mp;
            E = M + e;
            if (M > pi)
                E = M - e;
            elseif (-pi < M)
                if (M < 0)
                    E = M - e;
                end
            end
            %         //cout << " Input meanAnomaly = " << meanAnomaly << endl;
            for k=1:1:10
                E  = M + e*sin(E);
                %//cout << " Mean Anomaly Solution " << E << endl;
            end
            % //      10 rounds of Newton's root finding method based on the above "seed".
            for k=1:1:10
                Eprime      = 1.0 - e*cos(E);
                E           = E + (M - E + e*sin(E))./Eprime;
            end
            %KeplerInv = E;
            
            cosE  = cos(E);
            sinE  = sin(E);
            %E     = atan2(sinE,cosE)
            kfac  = 1.0 - e*cosE;
            cosnu = (cosE - e)./kfac;
            sinnu = rootfac*sinE./kfac;
            %     nu    = atan2(sinnu, cosnu)
            %     while(nu <   0.0)
            %         nu = nu + twopi;
            %     end
            %     while(nu > twopi)
            %         nu = nu - twopi;
            %     end
            %     sinnu = sin(nu);
            %     cosnu = cos(nu);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Pvec    = [ cosO*cosom - sinO*sinom*cosI;...
            %             sinO*cosom + cosO*sinom*cosI;...
            %             sinom*sinI];
            % Qvec    = [-cosO*sinom - sinO*cosom*cosI;...
            %            -sinO*sinom + cosO*cosom*cosI;...
            %             cosom*sinI];
            % Wvec    = [ sinO*sinI; -cosO*sinI; cosI];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Position and Unit Vector along the Position
            Rpos  = [];
            Xpos  = a*kfac.*(cosnu*Pvec(1) + sinnu*Qvec(1));
            Ypos  = a*kfac.*(cosnu*Pvec(2) + sinnu*Qvec(2));
            Zpos  = a*kfac.*(cosnu*Pvec(3) + sinnu*Qvec(3));
            Rpos  = DU*[Xpos' Ypos' Zpos'];
            %Rmag  = sqrt(Xpos.^2 + Ypos.^2 + Zpos.^2);
            %Runit = [(Xpos./Rmag)' (Ypos./Rmag)' (Zpos./Rmag)'];
            rtpinv   = sqrt(1.0./p);                           % line 53
            vorbitx  = rtpinv.*(-sinnu*Pvec(1) + (e + cosnu)*Qvec(1));     % line 54
            vorbity  = rtpinv.*(-sinnu*Pvec(2) + (e + cosnu)*Qvec(2));     % line 55
            vorbitz  = rtpinv.*(-sinnu*Pvec(3) + (e + cosnu)*Qvec(3));     % line 56
            Rdot  = VU*[vorbitx' vorbity' vorbitz'];
            %B = Runit';
            %B = B .^2;
            %sum(B);
            %LineCheck    = 50;                % Rx
            [H_I, Jacob, Hess] = WorkOrder(this, 50);
            Hessian{1}  = Hess;
            Jacobian      = [Jacobian; Jacob];   % row 1
            
            %LineCheck    = 51;                % Ry
            [H_I, Jacob, Hess] = WorkOrder(this, 51);
            Hessian{2}  = Hess;
            Jacobian      = [Jacobian; Jacob];   % row 1
            
            %LineCheck    = 52;                % Rz
            [H_I, Jacob, Hess] = WorkOrder(this, 52);
            Hessian{3}  = Hess;
            Jacobian      = [Jacobian; Jacob];   % row 1
            
            %LineCheck    = 54;                % Vx
            [H_I, Jacob, Hess] = WorkOrder(this, 54);
            Hessian{4}  = Hess;
            Jacobian      = [Jacobian; Jacob];   % row 1
            
            %LineCheck    = 55;                % Vy
            [H_I, Jacob, Hess] = WorkOrder(this, 55);
            Hessian{5}  = Hess;
            Jacobian      = [Jacobian; Jacob];   % row 1
            
            %LineCheck    = 56;                % Vz            
            [H_I, Jacob, Hess] = WorkOrder(this, 56);
            Hessian{6}  = Hess;
            Jacobian      = [Jacobian; Jacob];   % row 1
        end