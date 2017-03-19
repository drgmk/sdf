import io
import os
from datetime import datetime

import numpy as np

from jinja2 import Template
from bokeh.resources import CDN,INLINE

from . import plotting
from . import db
from . import templates
from . import utils
from . import config as cfg


def www_all(results,update=False):
    """Generate www material."""

    print(" Web")

    # see whether index.html needs updating (unless update enforced)
    index_file = results[0].path+'/index.html'
    
    if os.path.exists(index_file):
        mtime = []
        for r in results:
            mtime.append( r.mtime )

        if os.path.getmtime(index_file) > np.max(mtime):
            if not update:
                print("   no update needed")
                return
        print("   updating")
    else:
        print("   generating")

    index(results,file=index_file)


def index(results,file='index.html'):
    """Make index.html - landing page with an SED and other info."""

    script,div = plotting.sed_components(results)

    template = Template(templates.index)
    bokeh_js = CDN.render_js()
    bokeh_css = CDN.render_css()
#    bokeh_js = INLINE.render_js()
#    bokeh_css = INLINE.render_css()
    html = template.render(bokeh_js=bokeh_js,
                           bokeh_css=bokeh_css,
                           css=templates.css,
                           plot_script=script,
                           plot_div=div,
                           phot_file=os.path.basename(results[0].rawphot),
                           main_id=results[0].obs_keywords['main_id'],
                           spty=results[0].obs_keywords['sp_type'],
                           ra=results[0].obs_keywords['raj2000'],
                           dec=results[0].obs_keywords['dej2000'],
                           plx=results[0].obs_keywords['plx_value'],
                           xids=db.sdb_xids(results[0].obs_keywords['id']),
                           best_fit=results[0].main_results_text(),
                           corner=os.path.basename(results[0].pmn_dir)+'/'+\
                                    os.path.basename(results[0].corner_plot),
                           creation_time=datetime.utcnow().strftime("%d/%m/%y %X")
                           )

    with io.open(file, mode='w', encoding='utf-8') as f:
        f.write(html)

