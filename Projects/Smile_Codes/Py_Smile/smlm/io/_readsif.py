import numpy as np
from datetime import datetime, timedelta


# ---------- 底层读取工具 ----------
_WS = b' \t\r\n'

def _read_line(f):
    """等效 fgetl:读到 \n 为止,返回不含换行的字符串;文件末尾返回 None。"""
    chars = bytearray()
    got = False
    while True:
        c = f.read(1)
        if not c:
            break
        got = True
        if c == b'\n':
            break
        chars += c
    if not got:
        return None
    return chars.decode('latin-1').rstrip('\r')


def _skip_bytes(f, n):
    """等效 skipBytes。"""
    f.read(n)


def _skip_lines(f, n):
    """等效 skipLines。"""
    for _ in range(n):
        if _read_line(f) is None:
            raise ValueError('Inconsistent image header.')


def _scan(f, count, dtype=float):
    """
    等效 fscanf('%f'/'%d', count):跳过前导空白,逐个读取空白分隔的数字。
    读完最后一个数字后,把紧跟的非数字字符退回(seek),保持位置与 MATLAB 一致。
    """
    if dtype is float:
        numeric = set(b'0123456789+-.eEdD')
    else:
        numeric = set(b'0123456789+-')
    results = []
    while len(results) < count:
        # 跳过前导空白
        while True:
            pos = f.tell()
            c = f.read(1)
            if not c:
                break
            if c[0] not in _WS:
                f.seek(pos)          # 退回非空白字符
                break
        # 读取一个数字 token
        token = bytearray()
        while True:
            pos = f.tell()
            c = f.read(1)
            if not c:
                break
            if c[0] in numeric:
                token += c
            else:
                f.seek(pos)          # 退回分隔符
                break
        if not token:
            break
        text = token.decode('latin-1').replace('d', 'e').replace('D', 'E')
        results.append(float(text) if dtype is float else int(text))
    return results


def _read_string(f):
    """等效 readString:先读长度 n,再 fgetl 吃掉本行剩余,最后读 n 个字符。"""
    n = _scan(f, 1, float)
    if not n or n[0] < 0:
        raise ValueError('Inconsistent string.')
    n = int(n[0])
    if _read_line(f) is None:        # 对应条件里的 fgetl
        raise ValueError('Inconsistent string.')
    return f.read(n).decode('latin-1')


# ---------- 主函数 ----------
def readsif(file):
    """读取 Andor SIF 文件 (version 4.21~4.30),返回包含头信息与图像数据的 dict。"""
    info = {}
    with open(file, 'rb') as f:
        if _read_line(f) != 'Andor Technology Multi-Channel File':
            raise ValueError('Not an Andor SIF image file.')

        _skip_bytes(f, 8)
        o = _scan(f, 6, float)
        info['date'] = datetime(1970, 1, 1) + timedelta(seconds=o[4])  # MATLAB datenum -> datetime
        info['temperature'] = o[5]

        _skip_bytes(f, 10)
        o = _scan(f, 5, float)
        info['exposureTime'] = o[1]
        info['cycleTime'] = o[2]
        info['accumulateCycles'] = o[4]
        info['accumulateCycleTime'] = o[3]

        _skip_bytes(f, 2)
        o = _scan(f, 5, float)
        info['stackCycleTime'] = o[0]
        info['pixelReadoutTime'] = o[1]
        info['gainDAC'] = o[4]                    # EM gain

        _scan(f, 5, float)                        # 丢弃
        _skip_bytes(f, 4)
        o = _scan(f, 15, float)
        info['VertShift'] = o[12]
        info['PreAmpGain'] = o[14]

        info['SerialNumber'] = _scan(f, 1, int)[0]
        _scan(f, 9, float)                        # 丢弃
        o = _scan(f, 5, int)
        subver = o[1]
        info['SoftwareVersion'] = [o[0], o[1], o[2], o[3], o[4]]

        _read_line(f)                             # fgetl
        info['detectorType'] = _read_line(f).rstrip()
        info['detectorSize'] = _scan(f, 2, int)

        info['fileName'] = _read_string(f)
        csize = _scan(f, 2, int)
        _skip_bytes(f, csize[1] + 1)              # chunk = fread(...,uint8)
        _skip_lines(f, 28)

        info['frameAxis'] = _read_string(f)
        info['dataType'] = _read_string(f)
        info['imageAxis'] = _read_string(f)

        o = _scan(f, 16, int)
        # 与 MATLAB 完全一致的排布
        imageArea = np.array([[o[1], o[4], o[6]],
                              [o[3], o[2], o[5]]])
        frameArea = np.array([[o[10], o[13]],
                              [o[12], o[11]]])
        frameBins = np.array([o[15], o[14]])
        info['imageArea'] = imageArea
        info['frameArea'] = frameArea
        info['frameBins'] = frameBins

        # s = (1 + diff(frameArea)) ./ frameBins ;  diff 取列方向 (row2 - row1)
        diff_fa = frameArea[1, :] - frameArea[0, :]
        s = ((1 + diff_fa) / frameBins).astype(int)      # (帧高, 帧宽), 如 [256, 256]

        # z = 1 + diff(imageArea(5:6));  列优先线性索引 5,6 -> o[6], o[5]
        z = int(1 + (o[5] - o[6]))                       # 帧数

        _read_string(f)                          # 丢弃
        _skip_lines(f, z)
        if subver > 25:
            _skip_lines(f, z)
        _read_line(f)

        count = int(np.prod(s) * z)
        data = np.fromfile(f, dtype='<f4', count=count)  # '*single' -> float32 小端
        info['imageData'] = data.reshape((int(s[0]), int(s[1]), z), order='F')  # 列优先

    return info