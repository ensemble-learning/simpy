# -*- mode: python -*-
a = Analysis(['e_2_g03conf.py'],
             pathex=['D:\\Nut1\\code\\simpy\\src\\lib'],
             hiddenimports=[],
             hookspath=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=1,
          name=os.path.join('build\\pyi.win32\\e_2_g03conf', 'e_2_g03conf.exe'),
          debug=False,
          strip=None,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=True,
               name=os.path.join('dist', 'e_2_g03conf'))
