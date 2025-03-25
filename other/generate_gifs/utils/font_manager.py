import matplotlib.font_manager as fm
from PIL import ImageFont


class FontManager:
    """Manages fonts for the GIF processor"""

    def __init__(self):
        """Initialize font manager and discover available fonts"""
        self.fonts = self._discover_fonts()

    def _discover_fonts(self):
        """Discover all available fonts in the system"""
        fonts = {}

        # Get all font paths from matplotlib font manager
        for font in fm.fontManager.ttflist:
            try:
                # Create a normalized name (lowercase, no spaces)
                name = font.name.lower().replace(' ', '-')
                # Store the font path
                fonts[name] = font.fname
            except Exception:
                continue

        return fonts

    def get_font(self, font_name='default', size=30):
        """
        Get a font by name and size
        
        Args:
            font_name (str): Name of the font
            size (int): Font size
            
        Returns:
            PIL.ImageFont: Font object
        """
        try:
            # Handle default font
            if font_name == 'default':
                return ImageFont.load_default()

            # Try exact match
            if font_name in self.fonts:
                return ImageFont.truetype(self.fonts[font_name], size)

            # Try fuzzy match (e.g., "arial-bold" matches "arial")
            fuzzy_matches = [path for name, path in self.fonts.items()
                             if font_name in name and name != 'default']
            if fuzzy_matches:
                return ImageFont.truetype(fuzzy_matches[0], size)

            raise ValueError(f"Font '{font_name}' not found")

        except Exception as e:
            print(f"Error loading font {font_name}: {str(e)}")
            return ImageFont.load_default()

    def list_fonts(self):
        """List all available fonts"""
        return ['default'] + sorted(list(self.fonts.keys()))

    def search_fonts(self, query):
        """
        Search for fonts containing the query string
        
        Args:
            query (str): Search term
            
        Returns:
            list: List of matching font names
        """
        query = query.lower()
        matches = [name for name in self.fonts.keys() if query in name]
        return sorted(matches)
