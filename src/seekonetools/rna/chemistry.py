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
        'B': (1, 1),
        'L': (1, 1)
    },
    'nolinker':{
        'shift': True,
        'pattern': 'A',
        'structure': 'B8X8B8X10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
        'B': (1, 1)
    },
    "P3CBGB":{
         'shift': False,
         'pattern': 'A',
         'structure': 'B17U12',
         'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt'),),
         'B': (0, 0)
    }
}

