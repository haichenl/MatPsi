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

        %% natom, 1 double
        function varargout = natom(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('natom', this.objectHandle, varargin{:});
        end
        
        %% nbasis, 1 double
        function varargout = nbasis(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('nbasis', this.objectHandle, varargin{:});
        end
        
        %% nelec, 1 double
        function varargout = nelec(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('nelec', this.objectHandle, varargin{:});
        end
        
        %% coord, (natom, 3) matrix
        function varargout = coord(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('coord', this.objectHandle, varargin{:});
        end
		
		%% overlap, (nbasis, nbasis) matrix
        function varargout = overlap(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('overlap', this.objectHandle, varargin{:});
        end
        
        %% kinetic, (nbasis, nbasis) matrix
        function varargout = kinetic(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('kinetic', this.objectHandle, varargin{:});
        end
        
        %% potential, (nbasis, nbasis) matrix
        function varargout = potential(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('potential', this.objectHandle, varargin{:});
        end
        
        %% separated potential, (nbasis, nbasis, natom) 3-d matrix
        function varargout = potential_sep(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('potential_sep', this.objectHandle, varargin{:});
        end
        
        %% tei_ijkl, 1 double
        function varargout = tei_ijkl(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('tei_ijkl', this.objectHandle, varargin{:});
        end

    end
end