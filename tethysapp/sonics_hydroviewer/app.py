from tethys_sdk.base import TethysAppBase, url_map_maker
from tethys_sdk.app_settings import CustomSetting, SpatialDatasetServiceSetting

class SonicsHydroviewer(TethysAppBase):
    """
    Tethys app class for SONICS Hydroviewer.
    """

    name = 'SONICS Hydroviewer'
    index = 'home'
    icon = 'sonics_hydroviewer/images/sonics_hydroviewer_logo.png'
    package = 'sonics_hydroviewer'
    root_url = 'sonics-hydroviewer'
    color = '#008B76'
    description = ''
    tags = '"Hydrology", "SONICS", "Hydroviewer", "Peru"'
    enable_feedback = False
    feedback_emails = []
    controller_modules = [ "controllers" ]

    def spatial_dataset_service_settings(self):
        """
		Spatial_dataset_service_settings method.
		"""
        return (
            SpatialDatasetServiceSetting(
                name='main_geoserver',
                description='spatial dataset service for app to use (https://tethys2.byu.edu/geoserver/rest/)',
                engine=SpatialDatasetServiceSetting.GEOSERVER,
                required=True,
            ),
        )

    # def url_maps(self):
    #     """
    #     Add controllers
    #     """
    #     UrlMap = url_map_maker(self.root_url)

    #     url_maps = (
    #         UrlMap(
    #             name='home',
    #             url='sonics-hydroviewer',
    #             controller='sonics_hydroviewer.controllers.home'
    #         ),
    #         UrlMap(
    #             name='get_hydrographs',
    #             url='get-hydrographs',
    #             controller='sonics_hydroviewer.controllers.get_hydrographs'
    #         ),
    #         UrlMap(
    #             name='get_simulated_discharge_csv',
    #             url='get-simulated-discharge-csv',
    #             controller='sonics_hydroviewer.controllers.get_simulated_discharge_csv'
    #         ),
    #         UrlMap(
    #             name='get-time-series',
    #             url='get-time-series',
    #             controller='sonics_hydroviewer.controllers.get_time_series'
    #         ),
    #         UrlMap(
    #             name='get_forecast_data_csv',
    #             url='get-forecast-data-csv',
    #             controller='sonics_hydroviewer.controllers.get_forecast_data_csv'
    #         ),
    #     )

    #     return url_maps

    def custom_settings(self):
        return (
            CustomSetting(
                name='folder',
                type=CustomSetting.TYPE_STRING,
                description="Floder where the SONICS Forecast are stored",
                required=True,
                default='/home/tethys/sonics',
            ),
        )