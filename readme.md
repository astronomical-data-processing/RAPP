update:
================================
    Add english version of readme.md.
    A little bit of code modification
********************************
# Warning!
    Do not use IDE to debugging the program
    Because RAPP uses the multiprocessing module
********************************
# Simple photometry tutorial:
>##      First import module:

        from module import core

>##      Then create photometric object:

        app = core.APpipeline(targ='xxx',
                              expo_key='EXPOSE',
                              date_key='DATE',)

        There are three required parameters, the fits folder path, the key of the exposure, and the key of the date.
        These two keys need to be found in the header of fits file.
        After creating the app object, you can use it to get the following parameters:
            app.targ:   list    The list of fits path.
            app.mask:   ndarr   Mask image. Default is an incut ellipse.
            app.bias:   ndarr   Bias image. Default is 0.
            app.dark:   ndarr   Dark image. Default is 0.
            app.flat:   ndarr   Flat image. Default is 1.
        The suggestion here is to fill in all the parameters and the full version looks like this:
            app = ap.APpipeline(targ=,      required
                                bias=,
                                dark=,
                                flat=,
                                mask=,
                                expo_key=,  required
                                date_key=,  required
                                count=,
                                N=)

>##      Once the app is created, then do the information initialize:

        app.info_init()

        This step will create app.info
        app.info is a pandas.DataFrame
        The structure is:(jd, (geometric radius, center of mass), star)
    
>##      After information initialized, do the match:

        app.match()

        The match() will add app.shifts to the app object.
        app.shifts is a ndarray, datatype is complex number
        Each complex number represents the shift of each image to the reference image
    
>##      After matched, do aperture photometry:

        app.ap()

        The ap() will add app.table to the app object.如此便会给app产生一个新的参数app.table
        app.table is a dictionary
        The structure is:(stars, jd, (magnitude, error))
        There are 2 important parameters in ap():
            a:      tuple   Default is (1.2, 2.4, 3.6)
                (aperture ratio, background inner ratio, background outer ratio)
            gain:   float   Default is 1.

>##      After aperture photometry, we can save the informathin we want:

        app.save(result='folder')
        app.draw(result='folder')

        save() is used to save the csv file and graph it.
        draw() is used to draw the reference image and mark stars.
********************************
# Overlay photometry tutorial:
>##      First:

        from module import core

>##      Then creat app:

        app = core.APpipeline(data='xxx',
                              bias='xxx',
                              dark='xxx',
                              flat='xxx',
                              expo_key='EXPOSE',
                              date_key='DATE',
                              N=1,)

        Because here is overlay photometry, count USES the default value of 6, which doesn't need to be too high, for reasons explained later

>##      Then information initialize::

        app.info_init()
    
>##      Then match:
        app.match()
    
>##      Then here is the key point, do image combine:

        img = app.img_combine()

>##      Find stars in image:

        info = app.find_star(raw=img,
                             ref=True,
                             count=10)

        Here's why we didn't have to provide count when we created the app. Because the weak stars in the combined image are definitely more obvious than the original data, you can find more stars in the combined image.
        And also we get a new info from the combined image.

>##      Then aperture photometry, but it also need to throw info to ap() as an argument:

        app.ap(info)

>##      In the later process, the difference from the simple version is on draw():

        app.draw(result='folder',
                 ref=img,
                 info0=info)
        
        Then, the reference image will draw by combine image

>##     Last:

        app.save(result='folder')
