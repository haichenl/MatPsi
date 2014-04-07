% MatPsi: An interface between Matlab and Psi4 
classdef MatPsi < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance  
        function this = MatPsi(varargin)
            this.objectHandle = MatPsi_mex('new', varargin{:});
        end
		
		%% Copy Constructor 
		function this = MatPsiCopy(this, varargin)
            this.objectHandle = MatPsi_mex('MatPsiCopy', this.objectHandle, varargin{:});
        end
        
        %% Destructor - Destroy the C++ class instance 
        function delete(this)
            MatPsi_mex('delete', this.objectHandle);
        end

        %% Molecule properties 
        % natom, 1 double 
        function varargout = natom(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('natom', this.objectHandle, varargin{:});
        end
        
        % nelec, 1 double 
        function varargout = nelec(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('nelec', this.objectHandle, varargin{:});
        end
        
        % coord, (natom, 3) matrix 
        function varargout = coord(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('coord', this.objectHandle, varargin{:});
        end
        
        % Zlist, (natom, 1) matrix 
        function varargout = Zlist(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('Zlist', this.objectHandle, varargin{:});
        end
        
        % Enuc, 1 double 
        function varargout = Enuc(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('Enuc', this.objectHandle, varargin{:});
        end
        
        %% Basis set properties 
        % nbasis, 1 double 
        function varargout = nbasis(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('nbasis', this.objectHandle, varargin{:});
        end
        
        % func2center, (nbasis, 1) matrix 
        function varargout = func2center(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('func2center', this.objectHandle, varargin{:});
        end
        
        % func2am, (nbasis, 1) matrix 
        function varargout = func2am(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('func2am', this.objectHandle, varargin{:});
        end
		
        %% One-electron integrals 
		% overlap, (nbasis, nbasis) matrix 
        function varargout = overlap(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('overlap', this.objectHandle, varargin{:});
        end
        
        % kinetic, (nbasis, nbasis) matrix 
        function varargout = kinetic(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('kinetic', this.objectHandle, varargin{:});
        end
        
        % potential, (nbasis, nbasis) matrix 
        function varargout = potential(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('potential', this.objectHandle, varargin{:});
        end
        
        % atom-separated potential, (nbasis, nbasis, natom) 3-d matrix 
        function varargout = potential_sep(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('potential_sep', this.objectHandle, varargin{:});
        end
        
        % environment potential for a point charge, (nbasis, nbasis) matrix 
        function varargout = potential_zxyz(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('potential_zxyz', this.objectHandle, varargin{:});
        end
        
        %% Two-electron integrals 
        % tei_ijkl, 1 double 
        function varargout = tei_ijkl(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('tei_ijkl', this.objectHandle, varargin{:});
        end
        
        % tei_alluniq, (nuniq, 1) vector 
        function varargout = tei_alluniq(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('tei_alluniq', this.objectHandle, varargin{:});
        end
        
        % tei_alluniqJK, (nuniq, 1) vector 
        function varargout = tei_alluniqJK(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('tei_alluniqJK', this.objectHandle, varargin{:});
        end
        
        %% SCF related 
        % HFnosymmMO2G, (nbasis, nbasis) matrix 
        function varargout = HFnosymmMO2G(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('HFnosymmMO2G', this.objectHandle, varargin{:});
        end

    end
end