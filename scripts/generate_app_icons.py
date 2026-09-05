"""
Generate Android App Icons from high-resolution master art for ProteinScope
"""
import os
import math
from PIL import Image, ImageDraw, ImageFilter

ANDROID_RES_DIR = os.path.abspath("android/app/src/main/res")
MOBILE_ICONS_DIR = os.path.abspath("mobile_app/icons")

SIZES = {
    "mipmap-mdpi": 48,
    "mipmap-hdpi": 72,
    "mipmap-xhdpi": 96,
    "mipmap-xxhdpi": 144,
    "mipmap-xxxhdpi": 192,
}

FOREGROUND_SIZES = {
    "mipmap-mdpi": 108,
    "mipmap-hdpi": 162,
    "mipmap-xhdpi": 216,
    "mipmap-xxhdpi": 324,
    "mipmap-xxxhdpi": 432,
}

def draw_master_icon(size=1024, is_round=False, is_foreground=False):
    # 2x supersampling for ultra-crisp antialiasing
    canvas_size = size * 2
    img = Image.new("RGBA", (canvas_size, canvas_size), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    cx, cy = canvas_size // 2, canvas_size // 2
    r_outer = canvas_size * 0.44

    if not is_foreground:
        # Draw background squircle or circle
        bg = Image.new("RGBA", (canvas_size, canvas_size), (7, 11, 22, 255))
        bg_draw = ImageDraw.Draw(bg)

        # Subtle radial glow from center
        for r in range(int(r_outer * 1.2), 0, -12):
            alpha = int(35 * (1.0 - r / (r_outer * 1.2)))
            bg_draw.ellipse([cx - r, cy - r, cx + r, cy + r], fill=(14, 165, 233, alpha))

        # Mask background
        mask = Image.new("L", (canvas_size, canvas_size), 0)
        mask_draw = ImageDraw.Draw(mask)
        if is_round:
            mask_draw.ellipse([cx - r_outer, cy - r_outer, cx + r_outer, cy + r_outer], fill=255)
        else:
            rx = int(canvas_size * 0.22)
            mask_draw.rounded_rectangle([cx - r_outer, cy - r_outer, cx + r_outer, cy + r_outer], radius=rx, fill=255)

        img.paste(bg, (0, 0), mask)

        # Outer border
        border_color = (56, 189, 248, 60)
        if is_round:
            draw.ellipse([cx - r_outer, cy - r_outer, cx + r_outer, cy + r_outer], outline=border_color, width=6)
        else:
            rx = int(canvas_size * 0.22)
            draw.rounded_rectangle([cx - r_outer, cy - r_outer, cx + r_outer, cy + r_outer], radius=rx, outline=border_color, width=6)

    # Scale factor for logo elements
    scale = canvas_size / 512.0

    # 1. Orbital Reticle Rings
    draw.ellipse([cx - 185*scale, cy - 185*scale, cx + 185*scale, cy + 185*scale], outline=(56, 189, 248, 70), width=int(3*scale))
    draw.ellipse([cx - 140*scale, cy - 140*scale, cx + 140*scale, cy + 140*scale], outline=(99, 102, 241, 60), width=int(2*scale))

    # 2. Glowing Center Aura
    for r in range(int(110*scale), 0, -8):
        alpha = int(45 * (1.0 - r / (110*scale)))
        draw.ellipse([cx - r, cy - r, cx + r, cy + r], fill=(56, 189, 248, alpha))

    # 3. Double Helix Ribbon Knot Points
    # Compute parametric smooth lemniscate / ribbon
    points_violet = []
    points_cyan = []
    num_pts = 180

    for i in range(num_pts):
        t = (i / num_pts) * 2 * math.pi
        # Lemniscate of Bernoulli with 3D tilt
        denom = 1 + math.sin(t)**2
        x = cx + (125 * scale * math.cos(t)) / denom
        y = cy + (125 * scale * math.sin(t) * math.cos(t)) / denom
        # Offset for ribbon depth
        z = math.sin(t)

        if z > -0.1:
            points_cyan.append((x, y))
        if z < 0.1:
            points_violet.append((x, y))

    # Draw Ribbon Back (Violet/Indigo)
    stroke_w = int(24 * scale)
    for i in range(len(points_violet) - 1):
        p1 = points_violet[i]
        p2 = points_violet[i+1]
        prog = i / max(1, len(points_violet))
        r = int(168 * (1 - prog) + 99 * prog)
        g = int(85 * (1 - prog) + 102 * prog)
        b = int(247 * (1 - prog) + 241 * prog)
        draw.line([p1, p2], fill=(r, g, b, 230), width=stroke_w)

    # Draw Hydrogen Bond Rungs
    rungs = [
        (cx - 50*scale, cy - 30*scale, cx - 25*scale, cy + 45*scale),
        (cx, cy - 65*scale, cx, cy + 65*scale),
        (cx + 50*scale, cy - 30*scale, cx + 25*scale, cy + 45*scale),
    ]
    for x1, y1, x2, y2 in rungs:
        draw.line([x1, y1, x2, y2], fill=(56, 189, 248, 160), width=int(4*scale))

    # Draw Ribbon Front (Electric Cyan / Cobalt)
    for i in range(len(points_cyan) - 1):
        p1 = points_cyan[i]
        p2 = points_cyan[i+1]
        prog = i / max(1, len(points_cyan))
        r = int(56 * (1 - prog) + 2 * prog)
        g = int(189 * (1 - prog) + 132 * prog)
        b = int(248 * (1 - prog) + 199 * prog)
        draw.line([p1, p2], fill=(r, g, b, 255), width=stroke_w)

    # Secondary Arc Rings
    draw.arc([cx - 95*scale, cy - 95*scale, cx + 95*scale, cy + 95*scale], start=30, end=190, fill=(56, 189, 248, 120), width=int(5*scale))
    draw.arc([cx - 95*scale, cy - 95*scale, cx + 95*scale, cy + 95*scale], start=210, end=370, fill=(168, 85, 247, 120), width=int(5*scale))

    # 4. Atomic Satellite Nodes
    nodes = [
        (cx - 110*scale, cy + 20*scale, 14*scale, (56, 189, 248)),
        (cx + 110*scale, cy - 20*scale, 14*scale, (168, 85, 247)),
        (cx - 60*scale, cy - 80*scale, 10*scale, (56, 189, 248)),
        (cx + 60*scale, cy + 80*scale, 10*scale, (99, 102, 241)),
    ]
    for nx, ny, nr, col in nodes:
        draw.ellipse([nx - nr, ny - nr, nx + nr, ny + nr], fill=col)
        draw.ellipse([nx - nr*0.4, ny - nr*0.4, nx + nr*0.4, ny + nr*0.4], fill=(255, 255, 255, 240))

    # 5. Glowing Central Active Site Core Energy Atom
    core_r = 18 * scale
    for cr in range(int(core_r * 2), int(core_r), -2):
        a = int(60 * (1.0 - (cr - core_r)/core_r))
        draw.ellipse([cx - cr, cy - cr, cx + cr, cy + cr], fill=(245, 158, 11, a))
    draw.ellipse([cx - core_r, cy - core_r, cx + core_r, cy + core_r], fill=(245, 158, 11, 255))
    draw.ellipse([cx - core_r*0.45, cy - core_r*0.45, cx + core_r*0.45, cy + core_r*0.45], fill=(255, 255, 255, 255))

    # Downsample using high-quality Lanczos filter
    final_img = img.resize((size, size), Image.Resampling.LANCZOS)
    return final_img

def main():
    print("Generating ProteinScope high-resolution app icons...")

    # 1. PWA icon (512x512 PNG)
    pwa_icon = draw_master_icon(512, is_round=False)
    pwa_path = os.path.join(MOBILE_ICONS_DIR, "icon.png")
    pwa_icon.save(pwa_path, "PNG")
    print(f"[OK] Saved PWA icon: {pwa_path}")

    # 2. Android Mipmap standard and round icons
    for folder, dim in SIZES.items():
        target_dir = os.path.join(ANDROID_RES_DIR, folder)
        os.makedirs(target_dir, exist_ok=True)

        # Standard icon
        std_img = draw_master_icon(dim, is_round=False)
        std_path = os.path.join(target_dir, "ic_launcher.png")
        std_img.save(std_path, "PNG")

        # Round icon
        rnd_img = draw_master_icon(dim, is_round=True)
        rnd_path = os.path.join(target_dir, "ic_launcher_round.png")
        rnd_img.save(rnd_path, "PNG")

        print(f"[OK] Generated {folder} ({dim}x{dim})")

    # 3. Android Adaptive Foreground icons
    for folder, dim in FOREGROUND_SIZES.items():
        target_dir = os.path.join(ANDROID_RES_DIR, folder)
        fg_img = draw_master_icon(dim, is_round=False, is_foreground=True)
        fg_path = os.path.join(target_dir, "ic_launcher_foreground.png")
        fg_img.save(fg_path, "PNG")
        print(f"[OK] Generated Foreground {folder} ({dim}x{dim})")

    print("\nAll Android icons generated successfully!")

if __name__ == "__main__":
    main()
