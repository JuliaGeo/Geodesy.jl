"""
    geoid(geoid::String; goids_folder::String = "~/Documents/geoids")
returns geoid::GeoArray{} height file [relative to WGS 84 (EPSG::4979) ellipsoid]

If geoid file not found in `goids_folder` it will be downloaded from:
https://www.agisoft.com/downloads/geoids/ 

Models Converted from USA NGA data by agisoft under Public Domain license.
"""
function geoid(geoid::String; goids_folder::String="~/Documents/geoids")

    if geoid == "egm84"
        # EGM84 30' geoid model
        # WGS 84 (EPSG::4979) to EGM84 height (EPSG::5798)
        geoid_fn = "us_nga_egm84_30.tif"
    elseif geoid == "egm96"
        # EGM96 15' geoid model
        # WGS 84 (EPSG::4979) to EGM96 height (EPSG::5773
        geoid_fn = "us_nga_egm96_15.tif"
    elseif geoid == "egm2008"
        # EGM2008 1' geoid model
        # WGS 84 (EPSG::4979) to EGM2008 height (EPSG::3855)
        geoid_fn = "us_nga_egm2008_1.tif"
    else
        error("geoid not recognized, valid geoid names are \"egm84\", \"egm96\" and \"egm2008\"")
    end

    path2goid = joinpath(goids_folder, geoid_fn)

    # download if file does not exist
    if !isfile(path2goid)
        if !isdir(goids_folder)
            error("goids folder does not exist: $goids_folder")
        else
            url = joinpath("https://s3-eu-west-1.amazonaws.com/download.agisoft.com/gtg/", geoid_fn)
            printstyled("local copy of $geoid file not found, downloading from: $url \n"; color = :blue, bold = true)
            HTTP.download(url, path2goid)
        end
    end

    # not sure if this Type decleration helps at all, feel free to delete
    return GeoArrays.read(path2goid; masked=false)::GeoArray{Float32, ArchGDAL.RasterDataset{Float32, ArchGDAL.IDataset}} 
end