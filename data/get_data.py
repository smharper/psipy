from __future__ import print_function
import sys

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen


if __name__ == '__main__':
    block_size = 16384
    base_url = 'http://web.mit.edu/smharper/Public/psipy-data/'
    files = ['fine_azim_coarse_E.h5', 'fine_E_coarse_azim.h5']

    for fname in files:
        url = base_url + fname
        response = urlopen(url)

        if sys.version_info[0] < 3:
            file_size = int(response.info().getheaders('Content-Length')[0])
        else:
            file_size = response.length

        print('Downloading {0}... '.format(fname), end='')
        downloaded = 0
        with open(fname, 'wb') as fh:
            chunk = response.read(block_size)
            while chunk:
                fh.write(chunk)
                downloaded += len(chunk)
                status = '{0:10}  [{1:3.2f}%]'.format(downloaded,
                                                  downloaded * 100. / file_size)
                print(status + chr(8)*len(status), end='')
                chunk = response.read(block_size)
        print('')
