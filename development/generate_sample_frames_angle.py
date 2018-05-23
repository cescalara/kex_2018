import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors
from IPython.display import display
import os, errno

class SampleGenerator():
    """
    generate labelled simulated data of background and 
    a simple track model to use as input to ML classificaiton
    """

    def __init__(self, track_model = None):

        # constants
        
        # frame size
        self._n_row = 48
        self._n_col = 48
        
        # mean background level
        self._mu_bg = 1 # photon/pixel/GTU
        
        # number of frames to be generated
        self._n_frame = 100
        
        # number of angle categories
        self._n_angles = 4

        # initialise samples
        self._bg_frames = []
        self._track_frames = []
        
        # initialise labels
        self._bg_labels = []
        self._track_labels = []

        # track_model status
        self.track_model = track_model
        
        # norm
        self.norm = colors.Normalize(0, track_model.mu_c + track_model.sigma_c)
        
    def __enter__(self):
        return self

    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.track_model = None

        
    def _factors(self, n):
        return set(
            factor for i in range(1, int(n**0.5) + 1) if n % i == 0
            for factor in (i, n//i))


    def _display(self, frames, number, variant):
        """
        display a random sample of the generated frames 
        """

        print ("displaying a random sample of", number, "frames: ")

        # divide into subplots for nice display
        set_fac = SampleGenerator._factors(self, number)
        
        n_subplts_x = int(sorted(list(set_fac))[1])
        n_subplts_y = int(sorted(list(set_fac))[-2])

        # display
        fig, axarr = plt.subplots(n_subplts_x, n_subplts_y,
                                  figsize = (n_subplts_y, n_subplts_x))
        
        for x in range(n_subplts_x):
            for y in range(n_subplts_y):

                # pick a random frame number
                ran_frame = int(round((np.random.uniform(0, self._n_frame - 1, 1)).item()))
                axarr[x, y].imshow(frames[ran_frame], origin = "lower", cmap = "viridis", norm=self.norm)
                axarr[x,y].axis("off")
        
        # remove whitespace
        fig.subplots_adjust(hspace = 0.05, wspace = 0.05)
        fig.colorbar(axarr[0,0].imshow(frames[0], norm=self.norm), ax=axarr, orientation='vertical', fraction=.1)
        
        # save the plot as eps
        filename = variant+'.eps'
        plt.savefig(filename, format='eps', dpi=400)

        # show the plot
        display(fig)
        plt.close("all")

            
    def background(self, n_frame_in, n_angles_in):
        """
        generate _n_frame of background 
        """
        
        self._n_frame = n_frame_in
        self._n_angles = n_angles_in

        # sample from poisson to fill background frames
        samples = np.random.poisson(self._mu_bg, self._n_row * self._n_col * self._n_frame)
        self._bg_frames = np.reshape(samples, (self._n_frame, self._n_row, self._n_col)).tolist()

        # display some information regarding the generated frames
        print ("genarated", self._n_frame, "frames of background")
        
        # generate labels
        A = np.zeros((self._n_frame,), dtype=np.float64)
        A = np.reshape(A, (-1,1))
        
        for i in range(self._n_angles-1):
            addthis = np.zeros((self._n_frame,), dtype=np.float64)
            addthis = np.reshape(addthis, (-1,1))
            A = np.concatenate([A, addthis], axis=1)
            
        theones = np.ones((self._n_frame,), dtype=np.float64)
        theones = np.reshape(theones, (-1,1))
        A = np.concatenate([A, theones], axis=1)
        
        self._bg_labels = A
        
        # display some information regarding the generated labels
        print ("genarated", self._n_frame, "background labels")
        
        SampleGenerator._display(self, self._bg_frames, 10, 'bg')

        
    def tracks(self, n_frame_in, n_angles_in):
        """
        generate _n_frame of tracks based on the TrackModel class
        """
        from skimage.draw import line_aa
        
        self._n_frame = n_frame_in
        self._n_angles = n_angles_in

        # define the default track model, if not passed on __init__
        if self.track_model is None:
            self.track_model = TrackModel()
        
        A = np.zeros((self._n_frame,), dtype=np.float64)
        A = np.reshape(A, (-1,1))
        
        for i in range(self._n_angles):
            zerovec = np.zeros((self._n_frame,), dtype=np.float64)
            zerovec = np.reshape(zerovec, (-1,1))
            A = np.concatenate([A, zerovec], axis=1)
            
        counter = 0

        for frame in range(self._n_frame):
            track_frame = np.zeros((48, 48), dtype=np.uint8)

            # generate a track from the model
            self.track_model.generate_track()
            
            for w in range(1): 
                rr, cc, val = line_aa(self.track_model.start_position[0],
                                      self.track_model.start_position[1],
                                      self.track_model.end_position[0],
                                      self.track_model.end_position[1])
                counts_tmp = self.track_model.counts
            
                # give decreasing brightness
                for i in range(len(val)):
                    track_frame[rr[i], cc[i]] = val[i] * counts_tmp
                    if i%6 == 0:
                        counts_tmp = counts_tmp - 1
                    if counts_tmp < 0:
                        counts_tmp = 0

                        
            # add background
            samples = np.random.poisson(1, self._n_row * self._n_col)
            bg_frame = np.reshape(samples, (self._n_row, self._n_col))
            track_frame = track_frame + bg_frame
            self._track_frames.append(track_frame)
            
            
         # generate labels
            thetatemp = self.track_model.theta
            intervalSize = np.pi/self._n_angles
            
            for j in range(self._n_angles):
                if (thetatemp > j*intervalSize and thetatemp < (j+1)*intervalSize) or (thetatemp > j*intervalSize+np.pi and thetatemp < (j+1)*intervalSize+np.pi):
                    A[counter][j] = 1
                    
            self._track_labels = A  
            counter = counter + 1
                
        # print some information regarding the generated frames
        print ("generated", self._n_frame, "frames of tracks")
        
        # display some information regarding the generated labels
        print ("genarated", self._n_frame, "track labels")

        SampleGenerator._display(self, self._track_frames, 10, 'track')

        
    def save(self):
        """
        save the generated frames to a text file
        """
        import pickle
        
        
        if not os.path.exists("samples"):
            os.makedirs("samples")

        if self._bg_frames != []:
            with open("samples/sample_generator_bg.dat", "wb") as f:
                pickle.dump(self._bg_frames, f)
                
            with open("samples/sample_generator_bg_labels.dat", "wb") as f:
                pickle.dump(self._bg_labels, f)
        else:
            print ("Error: have not generated any background frames, so nothing to save")

        if self._track_frames != []:
            with open("samples/sample_generator_track.dat", "wb") as f:
                pickle.dump(self._track_frames, f)
                
            with open("samples/sample_generator_track_labels.dat", "wb") as f:
                pickle.dump(self._track_labels, f)
        else:
            print ("Error: have not generated any track frames, so nothing to save")

        print ("saved generated frames to file in samples/")
        


        
            
class TrackModel():
    """
    define a simple toy model to generate 
    UHECR - like tracks
    
    Parameters describing the track:

    start_position: pixel coordinates of starting position of the track
    length: in pixels
    width: in pixels
    theta: rotation angle (0 - 360 deg)
    phi: incidence angle (0, 45, 90 deg)
    counts: maximum # of counts in a single pixel
    """

    def __init__(self):

        #configurable parameters
        # start position pixel coordinates
        self.start_pos_min = 5
        self.start_pos_max = 43
        
        # length distribution (Gaussian)
        self.mu_l = 20
        self.sigma_l = 4
        
        # width distribution (Gaussian)
        self.mu_w = 3
        self.sigma_w = 1
        
        # angular distribtution (Uniform)
        self.theta_min = 0
        self.theta_max = 2 * np.pi
        self.phi_min = 0
        self.phi_max = np.pi / 2
        
        # counts distibution (Gaussian)
        self.mu_c = 15
        self.sigma_c = 3


        # track parameters
        self.start_position = np.zeros(2)
        self.end_position = np.zeros(2)
        self.length = 0
        self.width = 0
        self.theta = 0
        self.phi = 0
        self.counts = 0

        
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        """
    
    def generate_track(self):
        """
        build the track based on set parameters
        """
        
        # sample parameters from their distibutions
        self.start_position = np.around(np.random.uniform(self.start_pos_min, self.start_pos_max, 2)).astype(int)
        self.length = round((np.random.normal(self.mu_l, self.sigma_l, 1)).item())
        self.width = int(round((np.random.normal(self.mu_w, self.sigma_w, 1)).item()))
        self.theta = np.random.uniform(self.theta_min, self.theta_max, 1)
        self.phi = np.random.uniform(self.phi_min, self.phi_max, 1)
        self.counts = round((np.random.normal(self.mu_c, self.sigma_c, 1)).item())

        # calculate end position
        self.end_position[0] = self.start_position[0] + self.length * np.cos(self.theta)
        if self.end_position[0] > 45:
            self.end_position[0] = 45
        if self.end_position[0] < 0:
            self.end_position[0] = 0
        self.end_position[1] = self.start_position[1] + self.length * np.sin(self.theta)
        if self.end_position[1] > 45:
            self.end_position[1] = 45 
        if self.end_position[1] < 0:
            self.end_position[1] = 0
        self.end_position = np.around(self.end_position).astype(int)
        
