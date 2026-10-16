-- Frame list for the noise / depth analysis (processing_version v2, all exposure times and gain settings)
-- Run on the server shell:
--   psql -U gwuteam -p 5433 gwu --csv -f 00_query_frames.sql -o /data3/cwkim/exptimecal/noise_v2_frames.csv
SELECT
    p.id AS processed_id, p.raw_frame_id, p.filename AS processed_filename, p.cat_filename,
    s.unified_filename AS raw_filename, s.obstime, s.mjd, s.unit_id, s.filter_id, s.tile_id,
    s.object_name, s.object_type, s.is_test_observation, s.exptime, s.gain, s.egain,
    s.binning_x, s.binning_y, s.ccdtemp, s.set_ccdtemp, s.pixscale_x, s.pixscale_y,
    s.airmass, s.moon_phase, s.moon_sep, s.cloud_fraction,
    p.master_bias_id, p.master_dark_id, p.master_flat_id, p.reduction_status, p.sanity_flag, p.is_production,
    p.skyval, p.skysig, p.seeing, p.saturate, p.aper_4, p.aper_5, p.zp_4, p.zp_5, p.ezp_4, p.ezp_5, p.ul5_4, p.ul5_5
FROM survey_processedscienceframe AS p
JOIN survey_scienceframe AS s ON s.id = p.raw_frame_id
WHERE p.processing_version = 'v2'
ORDER BY s.obstime;
