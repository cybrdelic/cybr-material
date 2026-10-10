"""Process high-water resident memory, normalized to KiB on every platform.

Windows uses native PeakWorkingSetSize via psutil; macOS getrusage is bytes,
Linux getrusage is KiB. Fail if high-water accounting is unavailable.
"""
import sys
from types import SimpleNamespace

RUSAGE_SELF = 0

def peak_rss_mib():
    if sys.platform == 'win32':
        import psutil
        info = psutil.Process().memory_info()
        return info.peak_wset / 1024**2
    import resource
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value / (1024**2 if sys.platform == 'darwin' else 1024)

def getrusage(who):
    if who != RUSAGE_SELF:
        raise NotImplementedError('Only process high-water RSS is supported')
    return SimpleNamespace(ru_maxrss=peak_rss_mib()*1024)

def accounting_receipt():
    return {'platform':sys.platform, 'peak_rss_mib':peak_rss_mib(),
            'source':'Windows PeakWorkingSetSize' if sys.platform=='win32' else 'native getrusage',
            'normalized_ru_maxrss_unit':'KiB', 'includes_prior_process_allocations':True}
