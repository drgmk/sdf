"""Jinja2 templates for web pages."""

index=('''<!DOCTYPE html>
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
            <div id="target_head">
            
            <div class="search">
            <form style="display:inline-block;" action = "/~grant/sdb/search/search.php" method = "GET">
            <input type = "text" name = "id"/>
            <input type = "submit" value = "Search" />
            </form>
            </div>

            <div class="target_info">
            <h1>{{ main_id }}</h1>
                [

                {% if ra %}
                    &alpha;={{ (ra/15.0)|round(1) }}h, &delta;={{ dec|round(1) }}&deg;
                {% endif %}
                
                {% if spty %}
                    | {{ spty }}
                {% endif %}
                
                {% if plx %}
                    | d={{ (1000.0/plx)|round(1) }} pc
                {% endif %}

                ]
            </div>
            
            {% if xids %}
                <div class="xids">
                {% for id in xids %}
                    {% if loop.index is greaterthan 1 %}
                        ,
                    {% endif %}
                    {{ id }}
                {% endfor %}
                </div>
            {% endif %}
            
            <div class="links">
            <table>
            <tr><td>
                {% if phot_file %}
                Input:
                <a href="{{ phot_file }}" title="Input photometry file">data</a>
                {% endif %}
            </td>
            <td>
                Posteriors:
                <a href="{{ par_dist }}" title="Fitting posteriors" target="_blank">fitting</a>
                |
                <a href="{{ derived_dist }}" title="Derived parameters" target="_blank">derived</a>
            </td>
            
            <td>
            
            External:
            <a href="http://simbad.u-strasbg.fr/simbad/sim-basic?submit=SIMBAD+search&Ident={{ main_id }}" target="_blank" title="Simbad">simbad</a>
            
            {% if ra %}
                | <a href="http://cassis.sirtf.com/atlas/cgi/radec.py?ra={{ ra }}&dec={{ dec }}&radius=20" target="_blank" title="Cornell IRS Spectra">cassis</a>
            
                | <a href="http://irsa.ipac.caltech.edu/applications/finderchart/#id=Hydra_finderchart_finder_chart&RequestClass=ServerRequest&DoSearch=true&subsize=0.083&thumbnail_size=medium&sources=DSS,SDSS,twomass,WISE,IRIS&overlay_catalog=true&catalog_by_radius=true&iras_radius=240&sdss_radius=5&twomass_radius=5&wise_radius=5&one_to_one=_none_&dss_bands=poss1_blue,poss1_red,poss2ukstu_blue,poss2ukstu_red,poss2ukstu_ir&SDSS_bands=u,g,r,i,z&twomass_bands=j,h,k&wise_bands=1,2,3,4&UserTargetWorldPt={{ ra }};{{ dec }};EQ_J2000&projectId=finderchart&searchName=finder_chart&shortDesc=Finder%20Chart&isBookmarkAble=true&isDrillDownRoot=true&isSearchResult=true" target="_blank" title="IRSA Finder Chart query">finder</a>
                
                | <a href="http://sha.ipac.caltech.edu/applications/Spitzer/SHA/#id=SearchByPosition&RequestClass=ServerRequest&DoSearch=true&SearchByPosition.field.radius=0.13888889000000001&UserTargetWorldPt={{ ra }};{{ dec }};EQ_J2000&SimpleTargetPanel.field.resolvedBy=nedthensimbad&MoreOptions.field.prodtype=aor,pbcd&shortDesc=Position&isBookmarkAble=true&isDrillDownRoot=true&isSearchResult=true" target="_blank" title="Spitzer Heritage Archive query">spitzer</a>
                
                | <a href="https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html?searchQuery={{ ra }},{{ dec }}" target="_blank" title="MAST Archive">mast</a>
                
            {% endif %}
            </td></tr>
            </table>
            </div>
            
            {% if best_fit %}
                <div class="best_fit">
                <div class="best_fit_title">Best fit:</div>
                {% for m in best_fit %}
                <div class="best_fit_line">{{ m }}</div>
                {% endfor %}
                </div>
            {% endif %}

        </div><!-- end header section -->
        {% endif %}

        <table>
        <tr><td>
        {{ plot_script|indent(8) }}
        {{ plot_div|indent(8) }}
        </td></tr>
        </table>
        
        {% if creation_time %}
        <div class="footer">
            <hr>
            <p>generated {{ creation_time }} UTC</p>
        </div>
        {% endif %}
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
        
        {% if creation_time %}
        <div class="footer">
            <hr>
            <p>generated {{ creation_time }} UTC</p>
        </div>
        {% endif %}
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
        <div id=="target_head">
            <div>
                <h1>{{ title }}</h1>
            </div>
        </div>
        {% endif %}

        {{ plot_div|indent(8) }}
        {{ plot_script|indent(8) }}
        
        {% if creation_time %}
        <div class="footer">
            <hr>
            <p>generated {{ creation_time }} UTC</p>
        </div>
        {% endif %}
    </body>
</html>
''')

css=('''<link rel="stylesheet" type="text/css" href="https://fonts.googleapis.com/css?family=Lato">
    
        <link rel="stylesheet" type="text/css" href="http://camd21.ast.cam.ac.uk/~grant/sdb/style.css">''')

datatable=('''<!DOCTYPE html>
<html lang="en">
    <head>
    <meta charset="utf-8"/>
    <meta content="text/html;charset=UTF-8" http-equiv="Content-type"/>
    <style>
    table.dataTable {width: auto !important; margin: 0 !important;}
    .dataTables_filter, .dataTables_paginate {float: left !important; margin-left:1em}
    </style>
    <link href="https://cdn.datatables.net/1.10.9/css/jquery.dataTables.css" rel="stylesheet" type="text/css"/>
    <script src="https://code.jquery.com/jquery-1.11.3.min.js">
    </script>
    <script src="https://cdn.datatables.net/1.10.9/js/jquery.dataTables.min.js">
    </script>
    </head>
    
    {{ css }}
    
    <body class="sample_table">
    
    <div id="target_head">
    
    <div class="search">
        <form style="display:inline-block;" action = "/~grant/sdb/search/search.php" method = "GET">
        <input type = "text" name = "id"/>
        <input type = "submit" value = "Search" />
        </form>
    </div>

    <div class="target_info">
        <h1>{{ name }}</h1>
        [ <a href="hr.html">HR diagram</a>
        | <a href="fnuvsr.html">Flux/radius</a>
        | <a href="{{ name }}.xml">votable</a>
        ]
    </div>
    
    </div>
    
    <script>
    $(document).ready(function() {
        $('#{{ name }}').dataTable({
            "order": [],
            "scrollX": true,
            "iDisplayLength": 20,
            "aLengthMenu": [[10, 20, 50, 100, 500, 1000, -1], [10, 20, 50, 100, 500, 1000, 'All']],
        "pagingType": "full_numbers"
        });
    } );
    </script>
    
    {{ table }}
    
    {% if creation_time %}
    <div class="footer">
        <hr>
        <p>generated {{ creation_time }} UTC</p>
    </div>
    {% endif %}
    
    </body>
</html>
''')
