import io
import os

import numpy as np

from jinja2 import Template
from bokeh.resources import CDN

from . import plotting
from . import db
from . import templates
from .utils import SdfError
from . import config as cfg


def www_all(results,update=False):
    """Generate www material."""

    print(" Web")

    # see whether index.html needs updating (unless update is enforced)
    index_file = results[0].path+'/index.html'
    
    if os.path.exists(index_file):
        pkltime = []
        for r in results:
            pkltime.append( os.path.getmtime(r.pickle) )

        if os.path.getmtime(index_file) > np.max(pkltime):
            if not update:
                print("   no update needed")
                return
        print("   updating web")
    else:
        print("   generating web")

    index(results,file=index_file)


def index(results,tab_order=None,file='index.html'):
    """Make index.html - landing page with an SED and other info."""

    script,div = plotting.sed_components(results,tab_order=tab_order)

    info = db.sdb_info(results[0].id)
    if info is not None:
        sdbid,main_id,xids,ra,dec,dist = info
    else:
        sdbid,main_id,xids,ra,dec,dist = None,None,None,None,None

    template = Template(templates.sed)
    bokeh_js = CDN.render_js()
    bokeh_css = CDN.render_css()
    html = template.render(bokeh_js=bokeh_js,
                           bokeh_css=bokeh_css,
                           css=templates.css,
                           plot_script=script,
                           plot_div=div,
                           sdbid=sdbid,main_id=main_id,
                           ra=ra/15.0,dec=dec,dist=dist
                           )

    with io.open(file, mode='w', encoding='utf-8') as f:
        f.write(html)

