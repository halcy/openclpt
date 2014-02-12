/**
 * BMP output code. Stateful and terrible, but easy to use.
 *
 * (c) L. Diener 2010
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

// This header file deals with BMP's
FILE* bmp_file;
int line_pos = 0;
int line_max = 0;

#pragma pack ( 1 )
typedef struct bmp_header {
	uint8_t file_type[ 2 ]; // "BM"
	uint32_t file_size;
	uint16_t reserved_1; // 0
	uint16_t reserved_2; // 0
	uint32_t pixel_offset; // 54
	uint32_t header_size; // 40
	uint32_t x_size;
	uint32_t y_size;
	uint16_t planes; // 1
	uint16_t bpp; // 24
	uint32_t compression; // 0
	uint32_t image_size; // 0
	uint32_t x_ppm; // 0
	uint32_t y_ppm; // 0
	uint32_t used_colors; // 0
	uint32_t important_colors; // 0
} bmp_header;

void bmp_init( const char* file_name, int x_size, int y_size ) {
	// Make us a header.
	bmp_header file_head;
	file_head.file_type [ 0 ] = 'B';
	file_head.file_type [ 1 ] = 'M';
	file_head.file_size = x_size * y_size * 3 + 54;
	file_head.reserved_1 = 0;
	file_head.reserved_2 = 0;
	file_head.pixel_offset = 54;
	file_head.header_size = 40;
	file_head.x_size = x_size;
	file_head.y_size = y_size;
	file_head.planes = 1;
	file_head.bpp = 24;
	file_head.compression = 0;
	file_head.image_size = 0;
	file_head.x_ppm = 0;
	file_head.y_ppm = 0;
	file_head.used_colors = 0;
	file_head.important_colors = 0;

	// Write file header.
	bmp_file = fopen( file_name, "wb" );
	fwrite( (char*)(&file_head), 1, sizeof(file_head), bmp_file );

	line_max = x_size * 3;
	line_pos = 0;
}

void bmp_pixel(int r, int g,  int b) {
	char tmp_r = (char)r;
	char tmp_g = (char)g;
	char tmp_b = (char)b;
	fwrite( &tmp_b, 1, 1, bmp_file );
	fwrite( &tmp_g, 1, 1, bmp_file );
	fwrite( &tmp_r, 1, 1, bmp_file );

	line_pos+=3;
	if( line_pos == line_max ) {
		while( line_pos < line_max + line_max % 4 ) {
			fwrite( &tmp_b, 1, 1, bmp_file );
			line_pos++;
		}
		line_pos = 0;
	}
}

void bmp_close() {
	fflush( bmp_file );
	fclose( bmp_file );
}
