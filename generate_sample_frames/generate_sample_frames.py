import numpy as np
from matplotlib import pyplot as plt

class SampleGenerator():
    """
    generate labelled simulated data of background and 
    a simple track model to use as input to ML classificaiton
    """

    def __init__(self):

        # constants
        # frame size
        self._n_row = 48
        self._n_col = 48
        
        # mean background level
        self._mu_bg = 1 # photon/pixel/GTU
        
        # number of frames to be generated
        self._n_frame = 10

        # initialise samples
        self._bg_frames = []
        self._track_frames = []
        
    def __enter__(self):
        return self

    def _display(self, frames, number):
        """
        display a random sample of the generated frames 
        """
        print "displaying a random sample of", number, "background frames: "

        # plot to check
        plt.imshow(frames[0], origin = "lower")
        plt.colorbar()
        

    def background():
        """
        generate _n_frame of background 
        """

        # sample from poisson to fill background frames
        samples = np.random.poisson(self._mu_bg, self._n_row * self._n_col * self._n_frame)
        self._bg_frames = np.reshape(samples, (self._n_frame, self._n_row, self._n_col))

        # display some information regarding the generated frames
        print "genarated", self._n_frame, "of background"

        SampleGenerator._display(self, self._bg_frames, 10)

    def tracks():
        """
        generate _n_frame of tracks based on the TrackModel class
        """
        from skimage.draw import line_aa

        # define the track model
        track_model = TrackModel()
        for frame in range(self._n_frames):
            track_frame = np.zeros((48, 48), dtype=np.uint8)

            for w in range(1): 
                rr, cc, val = line_aa(start_position[0] + w, start_position[1], end_position[0] + w, end_position[1])
                counts_tmp = counts
            
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

                
        # print some information regarding the generated frames
        print "generated", self._n_frame, "of tracks"

        SampleGenerator.display(self, self._track_frames, 10)
        

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

    
    def _generate_tracks(self):
        """
        build the track based on set parameters
        """
        
        # sample parameters from their distibutions
        self.start_position = np.around(np.random.uniform(start_pos_min, start_pos_max, 2)).astype(int)
        self.length = round(np.random.normal(mu_l, sigma_l, 1))
        self.width = int(round(np.random.normal(mu_w, sigma_w, 1)))
        self.theta = np.random.uniform(theta_min, theta_max, 1)
        self.phi = np.random.uniform(phi_min, phi_max, 1)
        self.counts = round(np.random.normal(mu_c, sigma_c, 1))

        # calculate end position
        self.end_position[0] = start_position[0] + length * np.cos(theta)
        if self.end_position[0] > 45:
            self.end_position[0] = 45
        if self.end_position[0] < 0:
            self.end_position[0] = 0
        self.end_position[1] = start_position[1] + length * np.sin(theta)
        if self.end_position[1] > 45:
            self.end_position[1] = 45 
        if self.end_position[1] < 0:
            self.end_position[1] = 0
        self.end_position = np.around(self.end_position).astype(int)
        
