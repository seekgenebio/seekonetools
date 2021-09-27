import os

__srcdir = os.path.dirname(os.path.abspath(__file__))

CHEMISTRY = {
    'SO01V3':{
        'shift': True,
        'pattern': 'A',
        'structure': 'B8L8B8L10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
        'linker': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'Linker1.txt'),
                   os.path.join(__srcdir, 'barcode', 'SO01V3', 'Linker2.txt'),),
        'misB': (1, 0),
        'misL': (1, 0)
    }
}

