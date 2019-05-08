import h5py
from enum import Enum
from batch_iv_analysis.cell_struct import cell


class iv_file_type(Enum):
    UNSUPPORTED = 0  # everything else
    GENERIC = 1  # meant to be able to just read csv/tsv files
    MCGEHEE = 2  # mcgroup file format
    SNAITH_LEGACY = 3  # old snaith file format
    SNAITH = 4  # new snaith labview format
    H5 = 5  # github.com/mutovis/control-software


def detect_filetype(fullpath):
    filetype = iv_file_type.UNSUPPORTED
    return filetype


def read_file(fullPath):
    """
    Load iv data from a file
    and return a list of cell_data data classes
    """

    logMessages = StringIO()
    fileName, fileExtension = os.path.splitext(fullPath)

    # what we'll return
    ret_list = []
    ret = Object()

    isMcFile = False  # true if this is a McGehee iv file format
    isSnaithLegacyFile = False  # true if this is a Snaith legacy iv file format
    isSnaithFile = False  # true if this is a Snaith iv file format
    isMyFile = False  # true if this is a custom solar sim iv file format
    isH5 = False  # true when this is an hdf5 file

    print("Processing:", fileName, file=logMessages)
    if h5py.is_hdf5(fullPath):  # hdf5 file processing route
        isH5 = True
        try:
            h5 = h5py.File(fullPath, 'r')
        except Exception as inst:
            print("Failed to open {:} with error:".format(fullPath))
            print(inst)
            return

        fr = h5.attrs['Format Revision']
        if hasattr(fr, 'decode'):
            h5rev = h5.attrs['Format Revision'].decode()
        else:
            h5rev = str(h5.attrs['Format Revision'])
        print(
            "Found HDF5 solar sim data format revision {:s} data file".format(h5rev))

        this_ver = version.parse(h5rev)

        for substrate_str in list(h5.keys()):
            substrate = h5['/'+substrate_str]
            for pixel_str in list(substrate.keys()):
                ret = Object()
                pixel = substrate[pixel_str]

                ret.substrate = substrate_str
                ret.pixel = pixel_str
                if 'Intensity [suns]' in h5.attrs:
                    ret.suns = h5.attrs['Intensity [suns]']
                else:
                    ret.sunsA = h5.attrs['Diode 1 intensity [suns]']
                    ret.sunsB = h5.attrs['Diode 2 intensity [suns]']
                    # TODO: use the correct diode intensity for specific pixels instead of averaging the two diodes
                    ret.suns = (ret.sunsA + ret.sunsB)/2

                if this_ver <= version.parse('1.1.0'):
                    ret.area = float(pixel.attrs['area']) / 1e4  # in m^2
                else:
                    ret.area = pixel.attrs['area']
                ret.vsTime = False

                if 'ssPmax' in pixel.attrs:
                    ret.ssPmax = pixel.attrs['ssPmax']

                if 'Voc' in pixel.attrs:
                    ret.Voc = pixel.attrs['Voc']

                if 'Isc' in pixel.attrs:
                    ret.Isc = pixel.attrs['Isc']

                if 'Vmpp' in pixel.attrs:
                    ret.Vmpp = pixel.attrs['Vmpp']

                if 'Impp' in pixel.attrs:
                    ret.Impp = pixel.attrs['Impp']

                # this is all the i-v data
                if 'all_measurements' in pixel:
                    iv_data = pixel['all_measurements']
                else:
                    print('WARNING: Could not find any data in {:}'.format(
                        fileName), file=logMessages)
                    return

                if 'V_oc dwell' in iv_data.attrs:
                    ret.ssVoc = iv_data[iv_data.attrs['V_oc dwell']]

                if 'I_sc dwell' in iv_data.attrs:
                    ret.ssIsc = iv_data[iv_data.attrs['I_sc dwell']]

                if 'MPPT' in iv_data.attrs:
                    ret.mppt = iv_data[iv_data.attrs['MPPT']]

                # now we pick out regions of interest from the big i-v data set
                if 'Snaith' in iv_data.attrs:
                    # I_sc --> V_oc sweep
                    snaith_region = iv_data[iv_data.attrs['Snaith']]
                    ret.VV = np.array([e[0] for e in snaith_region])
                    ret.II = np.array([e[1] for e in snaith_region])
                    ret.reverseSweep = False
                    ret_list.append(copy.deepcopy(ret))
                if 'Sweep' in iv_data.attrs:
                    # V_oc --> I_sc sweep
                    sweep_region = iv_data[iv_data.attrs['Sweep']]
                    if ret_list != []:
                        ret_list.append(copy.deepcopy(ret_list[-1]))
                    else:
                        ret_list.append(copy.deepcopy(ret))
                    ret_list[-1].VV = np.array([e[0] for e in sweep_region])
                    ret_list[-1].II = np.array([e[1] for e in sweep_region])
                    ret_list[-1].reverseSweep = True

    else:  # (legacy) non-h5py file
        ret.reverseSweep = False
        if fileExtension == '.csv':
            delimiter = ','
        elif fileExtension == '.tsv':
            delimiter = '\t'
        else:
            delimiter = None

        fp = open(fullPath, mode='r')
        fileBuffer = fp.read()
        fp.close()
        min_length = 400
        if len(fileBuffer) < min_length:
            print('Could not read' + fileName +
                  '. This file is less than {:} characters long.'.format(min_length), file=logMessages)
            return
        head = fileBuffer[0:10]
        tail = fileBuffer[-min_length:-1]

        # mcFile test:
        # the first line is not a comment
        if (not head.__contains__('#')) and (head.__contains__('/')) and (head.__contains__('\t')):
            nMcHeaderLines = 25  # number of header lines in mcgehee IV file format
            # the first 8 chars do not contain comment symbol and do contain / and a tab, it's safe to assume mcgehee iv file format
            isMcFile = True
            # comment out the first 25 rows here
            fileBuffer = '#'+fileBuffer
            fileBuffer = fileBuffer.replace('\n', '\n#', nMcHeaderLines-1)
        # snaithFile test:
        elif 'suns:\t' in tail:
            nSnaithFooterLines = 11  # number of footer lines in snaith IV file format
            isSnaithLegacyFile = True
            delimiter = '\t'
            if (fileExtension == '.liv1') or (fileExtension == '.div1'):
                ret.reverseSweep = True
            if (fileExtension == '.liv2') or (fileExtension == '.div2'):
                ret.reverseSweep = False
            fileBuffer = fileBuffer[::-1]  # reverse the buffer
            # comment out the footer lines
            fileBuffer = fileBuffer.replace('\n', '#\n', nSnaithFooterLines+1)
            fileBuffer = fileBuffer[::-1]  # un-reverse the buffer
            fileBuffer = fileBuffer[:-3]  # remove the last (extra) '\r\n#'
        elif 'i-v file' in head:
            isMyFile = True
        elif '(# suns):\t' in tail:
            isSnaithFile = True
            footerLines = 21
            delimiter = '\t'
            if (fileExtension == '.liv1') or (fileExtension == '.div1'):
                ret.reverseSweep = True
            if (fileExtension == '.liv2') or (fileExtension == '.div2'):
                ret.reverseSweep = False
            fileBuffer = '#' + fileBuffer
            fileBuffer = fileBuffer[::-1]  # reverse the buffer
            # comment out the footer lines
            fileBuffer = fileBuffer.replace('\n', '#\n', footerLines+1)
            fileBuffer = fileBuffer[::-1]  # un-reverse the buffer
            fileBuffer = fileBuffer[:-3]  # remove the last (extra) '\r\n#'

        else:
            raise ValueError("Unknown input file type!")

        splitBuffer = fileBuffer.splitlines(True)

        ret.substrate = '?'
        ret.pixel = '?'
        ret.suns = 1
        ret.area = 1 * 1e-4  # in m^2
        ret.vsTime = False  # this is not an i,v vs t data file
        # extract comments lines and search for area and intensity
        comments = []
        for line in splitBuffer:
            if line.startswith('#'):
                comments.append(line)
                if 'Area' in line:
                    numbersHere = [
                        float(s) for s in line.split() if ivAnalyzer.isNumber(s)]
                    if len(numbersHere) is 1:
                        ret.area = numbersHere[0] * 1e-4
                elif 'I&V vs t' in line:
                    if float(line.split(' ')[5]) == 1:
                        ret.vsTime = True
                elif 'Number of suns:' in line or '(# suns)' in line:
                    numbersHere = [
                        float(s) for s in line.split() if ivAnalyzer.isNumber(s)]
                    if len(numbersHere) is 1:
                        ret.suns = numbersHere[0]
                elif line.startswith('#Pixel:'):
                    splitted = line.split('\t')
                    ret.pixel = splitted[1].strip()
                elif line.startswith('#Position:'):
                    splitted = line.split('\t')
                    ret.substrate = splitted[1].upper().strip()

        # for converstion to current density[mA/cm^2]
        jScaleFactor = 1000 / (ret.area*1e4)

        c = StringIO(fileBuffer)  # makes string look like a file

        #read in data
        try:
            data = np.loadtxt(c, delimiter=delimiter)
        except:
            print('Could not read' + fileName +
                  '. Prepend # to all non-data lines and try again', file=logMessages)
            return
        if data.size == 0:
            print('WARNING: Could not find any data in {:}'.format(
                fileName), file=logMessages)
            return
        if isMyFile:
            ret.VV = data[:, 2]
            ret.II = data[:, 3]
        elif isSnaithFile:
            ret.VV = data[:, 1]
            ret.II = data[:, 2]
        else:
            ret.VV = data[:, 0]
            ret.II = data[:, 1]
        if isMcFile or isSnaithLegacyFile:  # convert from current density to amps through soucemeter
            ret.II = ret.II/jScaleFactor

        ret_list.append(ret)

    for i in range(len(ret_list)):
        # prune data points that share the same voltage
        u, indices = np.unique(ret_list[i].VV, return_index=True)
        ret_list[i].VV = ret_list[i].VV[indices]
        ret_list[i].II = ret_list[i].II[indices]

        # sort data by ascending voltage
        newOrder = ret_list[i].VV.argsort()
        ret_list[i].VV = ret_list[i].VV[newOrder]
        ret_list[i].II = ret_list[i].II[newOrder]

    return ret_list
