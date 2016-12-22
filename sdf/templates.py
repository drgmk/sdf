"""Jinja2 templates for web pages."""

sed=('''<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8">
        <title>{{ main_id }}</title>
        {{ bokeh_css }}
        {{ bokeh_js }}

        {{ css }}
        
    </head>
    
    <body class="sed">
    
        {% if main_id %}
        <table class="head_table"><tr><td>
            <h1>{{ main_id }}</h1>
            [ <a href="http://simbad.u-strasbg.fr/simbad/sim-basic?submit=SIMBAD+search&Ident={{ main_id }}" target="_blank" \>simbad</a>

            {% if ra %}
            | <a href="http://irsa.ipac.caltech.edu/applications/finderchart/#id=Hydra_finderchart_finder_chart&RequestClass=ServerRequest&DoSearch=true&subsize=0.083&thumbnail_size=medium&sources=DSS,SDSS,twomass,WISE,IRIS&overlay_catalog=true&catalog_by_radius=true&iras_radius=240&sdss_radius=5&twomass_radius=5&wise_radius=5&one_to_one=_none_&dss_bands=poss1_blue,poss1_red,poss2ukstu_blue,poss2ukstu_red,poss2ukstu_ir&SDSS_bands=u,g,r,i,z&twomass_bands=j,h,k&wise_bands=1,2,3,4&UserTargetWorldPt={{ ra }};{{ dec }};EQ_J2000&projectId=finderchart&searchName=finder_chart&shortDesc=Finder%20Chart&isBookmarkAble=true&isDrillDownRoot=true&isSearchResult=true" target="_blank">finder</a>
            {% endif %}
            
            {% if sdbid %}
            | {{ sdbid }}
            {% endif %}
            
            ]
        </td>
        
        <td class="right">
            <a href="/~grant/sdb/">home</a>
            <form style="display:inline-block;" action = "/~grant/sdb/search/search.php" method = "GET">
            <input type = "text" name = "id"/>
            <input type = "submit" value = "Search" />
            </form>
        </td></tr></table>
        
        {% if xids %}
        <div>
            {% for id in xids %}
                <span>{{ id }}</span>
            {% endfor %}
        </div>
        {% endif %}

        {% endif %}

        {{ plot_div|indent(8) }}
        {{ plot_script|indent(8) }}
        
    </body>
</html>
''')

generic=('''<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8">
        <title>{{ title }}</title>
        {{ bokeh_css }}
        {{ bokeh_js }}

        {{ css }}
        
    </head>
    
    <body>
    
        {% if title %}
        <h1>{{ title }}</h1>
        {% endif %}

        {{ plot_div|indent(8) }}
        {{ plot_script|indent(8) }}
        
    </body>
</html>
''')

generic_wide=('''<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8">
        <title>{{ title }}</title>
        {{ bokeh_css }}
        {{ bokeh_js }}

        {{ css }}
        
    </head>
    
    <body class="wide">
    
        {% if title %}
        <h1>{{ title }}</h1>
        {% endif %}

        {{ plot_div|indent(8) }}
        {{ plot_script|indent(8) }}
        
    </body>
</html>
''')

css=('''<link rel="stylesheet" type="text/css" href="https://fonts.googleapis.com/css?family=Lato">
    
        <style>
            html {
                width: 100%;
                height: 100%;
            }
            body {
                font-family: Lato;
                width: 850px;
                height: 100%;
                margin: auto;
            }
            body.wide {
                width: 100%
            }
            .head_table {
                width:100%;
            }
            table.head_table td.right {
                vertical-align: top;
                text-align: right;
            }
            h1 {
                display: inline;
                padding: 0px 10px 0px 0px;
            }
        </style>''')
