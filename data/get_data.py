import urllib2


if __name__ == '__main__':
    base_url = 'http://web.mit.edu/smharper/Public/psipy-data/'
    files = ['fine_azim_coarse_E.h5', 'fine_E_coarse_azim.h5']

    for fname in files:
        url = base_url + fname
        response = urllib2.urlopen(url)

        with open(fname, 'wb') as fh:
            fh.write(response.read())
