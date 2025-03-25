from PIL import Image, ImageDraw
from io import BytesIO
import os
from IPython.display import display, Image as IPImage

from .font_manager import FontManager


class GIFProcessor:
    """Processes GIF files with various transformations and overlays"""

    def __init__(self, gif_path):
        """
        Initialize with path to GIF file
        
        Args:
            gif_path (str): Path to the GIF file
        """
        self.gif = Image.open(gif_path)
        self.frames = []
        self.durations = []

        # Extract all frames and their durations
        try:
            while True:
                self.frames.append(self.gif.copy())
                self.durations.append(self.gif.info.get('duration', 100))
                self.gif.seek(self.gif.tell() + 1)
        except EOFError:
            pass

    def resize(self, width, height):
        """
        Resize all frames to specified dimensions
        
        Args:
            width (int): Target width
            height (int): Target height
            
        Returns:
            GIFProcessor: Self for method chaining
        """
        self.frames = [frame.resize((width, height), Image.Resampling.LANCZOS)
                       for frame in self.frames]
        return self

    def make_square(self, size=None):
        """
        Crop the GIF to a square from the center.
        If size is provided, the output will be resized to size x size.
        If size is None, the square will be sized to the smaller dimension.
        
        Args:
            size (int, optional): Size of the square output
            
        Returns:
            GIFProcessor: Self for method chaining
        """
        if not self.frames:
            return self

        # Get dimensions from first frame
        width, height = self.frames[0].size

        # Calculate crop box for square
        if width > height:
            # Landscape orientation
            left = (width - height) // 2
            top = 0
            right = left + height
            bottom = height
        else:
            # Portrait orientation
            left = 0
            top = (height - width) // 2
            right = width
            bottom = top + width

        # Apply crop to all frames
        self.frames = [frame.crop((left, top, right, bottom)) for frame in self.frames]

        # Resize if size is specified
        if size is not None:
            self.frames = [frame.resize((size, size), Image.Resampling.LANCZOS)
                           for frame in self.frames]

        return self

    def add_text(self, text, position, font_path=None, font_size=30,
                 color=(255, 255, 255), stroke_width=2, stroke_fill=(0, 0, 0)):
        """
        Add text overlay to all frames
        
        Args:
            text (str): Text to add
            position (tuple): Position (x, y) to place text
            font_path (str): Font name or path
            font_size (int): Font size
            color (tuple): RGB color for text
            stroke_width (int): Width of text outline
            stroke_fill (tuple): RGB color for text outline
            
        Returns:
            GIFProcessor: Self for method chaining
        """
        font_manager = FontManager()
        font = font_manager.get_font(font_name=font_path, size=font_size)

        for i, frame in enumerate(self.frames):
            # Convert to RGBA before drawing
            frame_rgba = frame.convert('RGBA')
            draw = ImageDraw.Draw(frame_rgba)
            draw.text(position, text, font=font, fill=color,
                      stroke_width=stroke_width, stroke_fill=stroke_fill)
            self.frames[i] = frame_rgba
        return self

    def add_image_overlay(self, overlay_path, position):
        """
        Add image overlay to all frames
        
        Args:
            overlay_path (str): Path to overlay image
            position (tuple): Position (x, y) to place overlay
            
        Returns:
            GIFProcessor: Self for method chaining
        """
        overlay = Image.open(overlay_path).convert('RGBA')

        for i, frame in enumerate(self.frames):
            frame_rgba = frame.convert('RGBA')
            frame_rgba.paste(overlay, position, overlay)
            self.frames[i] = frame_rgba
        return self

    def optimize(self, quality=100):
        """
        Optimize GIF by reducing colors
        
        Args:
            quality (int): Quality from 1 (worst) to 100 (best)
        """
        if not self.frames:
            return

        # Convert frames to RGB mode for saving
        rgb_frames = []
        for frame in self.frames:
            rgb_frame = frame.convert('RGB')
            # Calculate number of colors based on quality
            n_colors = max(min(256, int(256 * (quality / 100))), 2)

            # Convert to P mode (palette) with optimized palette
            rgb_frame = rgb_frame.quantize(
                colors=n_colors,
                method=Image.Quantize.MEDIANCUT,
                dither=Image.Dither.FLOYDSTEINBERG
            )
            rgb_frames.append(rgb_frame)

        self.frames = rgb_frames

    def save(self, output_path, optimize=False, quality=100):
        """
        Save the processed GIF with optimization options
        
        Args:
            output_path (str): Path to save the GIF
            optimize (bool): Whether to optimize the GIF
            quality (int): Quality from 1 (worst) to 100 (best).
                         Lower quality means smaller file size.
        """
        if not self.frames:
            return

        if optimize:
            self.optimize(quality)

        rgb_frames = [frame.convert('RGB') for frame in self.frames]

        # Save with optimization
        rgb_frames[0].save(
            output_path,
            save_all=True,
            append_images=rgb_frames[1:],
            optimize=optimize,
            duration=self.durations,
            loop=0,
            format='GIF',
            # Additional optimization parameters
            disposal=2,  # Clear the frame before rendering the next
            quality=quality
        )
        print("Size on disk:", f"{os.path.getsize(output_path) / 1024 / 1024:.2f} MB")

    def display(self):
        """Display the current state of the GIF in the notebook"""
        temp_path = '_temp_display.gif'
        self.save(temp_path)
        display(IPImage(filename=temp_path))
        os.remove(temp_path)

    def get_size(self, optimize=False, quality=100):
        """
        Get the size of the processed GIF in bytes without saving to disk
        
        Args:
            optimize (bool): Whether to optimize the GIF
            quality (int): Quality from 1 (worst) to 100 (best)
            
        Returns:
            int: Size in bytes
        """
        if self.frames:
            # Convert frames back to RGB mode for saving
            rgb_frames = [frame.convert('RGB') for frame in self.frames]
            with BytesIO() as buffer:
                rgb_frames[0].save(
                    buffer,
                    save_all=True,
                    append_images=rgb_frames[1:],
                    optimize=optimize,
                    quality=quality,
                    duration=self.durations,
                    loop=0,
                    format='GIF'
                )
                return buffer.tell()
        return 0
