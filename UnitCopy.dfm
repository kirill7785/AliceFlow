object FormCopyObject: TFormCopyObject
  Left = 336
  Top = 183
  Caption = 'Copy object'
  ClientHeight = 232
  ClientWidth = 226
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnClose = FormClose
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object PanelMain: TPanel
    Left = 8
    Top = 8
    Width = 217
    Height = 217
    Color = clMoneyGreen
    TabOrder = 0
    object Labelnc: TLabel
      Left = 16
      Top = 16
      Width = 83
      Height = 13
      Caption = 'Number of copies'
    end
    object Editnc: TEdit
      Left = 112
      Top = 8
      Width = 73
      Height = 21
      TabOrder = 0
      Text = '1'
    end
    object GroupBoxTranslate: TGroupBox
      Left = 16
      Top = 40
      Width = 185
      Height = 121
      Caption = 'Translate'
      Color = clSkyBlue
      ParentColor = False
      TabOrder = 1
      object Labelx: TLabel
        Left = 8
        Top = 24
        Width = 33
        Height = 13
        Caption = 'Xoffset'
      end
      object Labely: TLabel
        Left = 8
        Top = 56
        Width = 33
        Height = 13
        Caption = 'Yoffset'
      end
      object Labelz: TLabel
        Left = 8
        Top = 88
        Width = 33
        Height = 13
        Caption = 'Zoffset'
      end
      object EditX: TEdit
        Left = 56
        Top = 16
        Width = 113
        Height = 21
        TabOrder = 0
        Text = '0.0'
      end
      object EditY: TEdit
        Left = 56
        Top = 48
        Width = 113
        Height = 21
        TabOrder = 1
        Text = '0.0'
      end
      object EditZ: TEdit
        Left = 56
        Top = 80
        Width = 113
        Height = 21
        TabOrder = 2
        Text = '0.0'
      end
    end
    object ButtonApply: TButton
      Left = 112
      Top = 176
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 2
      OnClick = ButtonApplyClick
    end
  end
end
